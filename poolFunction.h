namespace pool {

	/*===== 0 =================== ARWHEAD Function ================ 5000 ===================*/
	double ARWHEADRanged(double *x, double *g, double *sum, INT index, INT begin, INT end) {
		double fx = 0.0;
		double group1(0.0);
		for (int i = begin; i < end - 1; i++)
		{
			group1 = x[i] * x[i] + x[index] * x[index];
			fx += group1 * group1 + 3.0 - 4.0*x[i];
			g[i] = -4.0 + 4.0*group1*x[i];
			*sum += 4.0*group1*x[index];
		}
		return fx;
	}
	double ARWHEAD
	(
		double *x,
		double *g,
		const int n
	)
	{
		g[n - 1] = 0;
		double fx = valgradPoolOneCongestion(ARWHEADRanged, x, g, n, n - 1);
		double group1(0.0);
		for (int i = block_size - 1; i < n - 1; i = i + block_size)
		{
			group1 = x[i] * x[i] + x[n - 1] * x[n - 1];
			fx += group1 * group1 + 3.0 - 4.0*x[i];
			g[i] = -4.0 + 4.0*group1*x[i];
			g[n - 1] += 4.0*group1*x[n - 1];
		}
		return fx;
	}
	void InitializeARWHEAD
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = 1.0;
	}
	/*===== 1 ================= BDQRTIC Function ============= 5000 ==================*/

	double BDQRTICRanged(double *x, double *g, double *sum, INT index, INT begin, INT end) {
		double fx(0.0);
		double group1(0.0), group2(0.0);
		for (int i = begin; i < end; i++)
			g[i] = 0;
		double current_squared[4] = { 0.,pow(x[0], 2),pow(x[1], 2),pow(x[2], 2) };
		double last(pow(x[index], 2));
		for (int i = begin; i < end - 4; i++)
		{
			current_squared[0] = current_squared[1];
			current_squared[1] = current_squared[2];
			current_squared[2] = current_squared[3];
			current_squared[3] = pow(x[i + 3], 2);
			group1 = 3.0 - 4.0*x[i];
			group2 = current_squared[0] + 2 * current_squared[1] + 3 * current_squared[2]
				+ 4 * current_squared[3] + 5 * last;
			fx += group1 * group1 + group2 * group2;
			g[i] += -8.0*group1 + 4.0*group2*x[i];
			g[i + 1] += 8.0*group2*x[i + 1];
			g[i + 2] += 12.0*group2*x[i + 2];
			g[i + 3] += 16.0*group2*x[i + 3];
			*sum += 20.0*group2*x[index];
		}
		return fx;
	}

	double BDQRTICHelper(double *x, double *g, int n, int i) {
		double group1 = 3.0 - 4.0*x[i];
		double group2 = pow(x[i], 2) + 2 * pow(x[i + 1], 2) + 3 * pow(x[i + 2], 2)
			+ 4 * pow(x[i + 3], 2) + 5 * pow(x[n - 1], 2);
		double fx = group1 * group1 + group2 * group2;
		g[i] += -8.0*group1 + 4.0*group2*x[i];
		g[i + 1] += 8.0*group2*x[i + 1];
		g[i + 2] += 12.0*group2*x[i + 2];
		g[i + 3] += 16.0*group2*x[i + 3];
		g[n - 1] += 20.0*group2*x[n - 1];
		return fx;
	}

	double BDQRTIC(double *x, double *g, const int n) {
		double fx = valgradPoolOneCongestion(BDQRTICRanged, x, g, n, n - 1);
		for (int i = block_size - 4; i < n - 4; i = i + block_size) {
			fx += BDQRTICHelper(x, g, n, i);
			fx += BDQRTICHelper(x, g, n, i + 1);
			fx += BDQRTICHelper(x, g, n, i + 2);
			fx += BDQRTICHelper(x, g, n, i + 3);
		}
		return fx;
	}

	void InitializeBDQRTIC(double *x, const int n) {
		for (int i = 0; i < n; i++)
			x[i] = 1.0;
	}


	/*===== 2 ================= BROYDN7D Function ============= 5000==================*/
	//precondition: n>=4
	double BROYDN7D
	(
		double *x,
		double *g,
		const int n
	)
	{
		double first(-2.0*x[1] + 1 + (3. - 2.0*x[0])*x[0]);
		double last(-x[n - 2] + 1 + (3. - 2.0*x[n - 1])*x[n - 1]);
		double fx = pow(fabs(first), 7 / 3.0) + pow(fabs(last), 7 / 3.0);
		double powFabsFirst4over3, powFabsLast4over3;
		for (int i = 0; i < n; i++)
			g[i] = 0;

		g[0] = (7.0 / 3)* pow(fabs(first), 4 / 3.0)*((first > 0) ? (1) : (-1))*(3. - 4.*x[0]);
		g[1] = -(14.0 / 3)* pow(fabs(first), 4 / 3.0)*((first > 0) ? (1) : (-1));
		g[n - 2] = -(7.0 / 3)* pow(fabs(last), 4 / 3.0)*((last > 0) ? (1) : (-1));
		g[n - 1] = (7.0 / 3)* pow(fabs(last), 4 / 3.0)*((last > 0) ? (1) : (-1))*(3. - 4.*x[n - 1]);

		last = x[0] + x[n / 2];
		fx += pow(fabs(last), 7 / 3.0);
		g[0] += (7.0 / 3)* pow(fabs(last), 4 / 3.0)*((last > 0) ? (1) : (-1));
		g[n / 2] += (7.0 / 3)* pow(fabs(last), 4 / 3.0)*((last > 0) ? (1) : (-1));

		for (int i = 1; i < n / 2; i++)
		{
			first = 1 - x[i - 1] - 2.0*x[i + 1] + (3. - 2.0*x[i])*x[i];
			last = x[i] + x[i + n / 2];
			fx += pow(fabs(first), 7 / 3.0) + pow(fabs(last), 7 / 3.0);

			powFabsFirst4over3 = (7.0 / 3)* pow(fabs(first), 4 / 3.0)* ((first > 0) ? (1) : (-1));
			powFabsLast4over3 = (7.0 / 3)* pow(fabs(last), 4 / 3.0) * ((last > 0) ? (1) : (-1));

			g[i - 1] += -powFabsFirst4over3;
			g[i] += powFabsFirst4over3 * (3 - 4 * x[i])
				+ powFabsLast4over3;
			g[i + 1] += -2 * powFabsFirst4over3;
			g[i + n / 2] += powFabsLast4over3;
		}
		for (int i = n / 2; i < n - 1; i++)
		{
			first = 1 - x[i - 1] - 2.0*x[i + 1] + (3. - 2.0*x[i])*x[i];
			fx += pow(fabs(first), 7 / 3.0);
			powFabsFirst4over3 = (7.0 / 3)* pow(fabs(first), 4 / 3.0)* ((first > 0) ? (1) : (-1));
			g[i - 1] += -powFabsFirst4over3;
			g[i] += powFabsFirst4over3 * (3 - 4 * x[i]);
			g[i + 1] += -2 * powFabsFirst4over3;
		}
		return fx;
	}
	void InitializeBROYDN7D
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = 1.0;
	}
	/*===== 3 =================== BRYBND Function ================ 5000 ===================*/
	//i < n - must be
	double BRYBND
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = 0.0;
		double group(0.0), next(x[1] * (1 + x[1])), secGroup(next), current(x[0] * (1 + x[0]));
		for (int i = 0; i < n; i++)
			g[i] = 0;
		group = x[0] * (2.0 + 5.0* x[0] * x[0]) + 1 - secGroup;
		fx += group * group;
		g[0] += 2.0*group*(2.0 + 15.0*x[0] * x[0]);
		g[1] -= 2.0*group*(1 + 2.0*x[1]);
		for (int i = 1; i < 6; i++)
		{
			secGroup -= next;
			secGroup += current;
			current = next;
			next = x[i + 1] * (1 + x[i + 1]);
			secGroup += next;
			group = x[i] * (2.0 + 5.0* x[i] * x[i]) + 1 - secGroup;
			fx += group * group;
			g[i] += 2.0*group*(2.0 + 15.0 * x[i] * x[i]);

			int  lo((0 > i - 5) ? (0) : (i - 5));
			for (int j = lo; j < i; j++)
				g[j] -= 2.0*group*(1 + 2.0*x[j]);
			g[i + 1] -= 2.0*group*(1 + 2.0*x[i + 1]);
		}
		for (int i = 6; i < n - 1; i++)
		{
			secGroup -= next;
			secGroup += current;
			current = next;
			next = x[i + 1] * (1 + x[i + 1]);
			secGroup += next;
			secGroup -= x[i - 6] * (1 + x[i - 6]);
			group = x[i] * (2.0 + 5.0* x[i] * x[i]) + 1 - secGroup;
			fx += group * group;
			g[i] += 2.0*group*(2.0 + 15.0 * x[i] * x[i]);

			for (int j = i - 5; j < i; j++)
				g[j] -= 2.0*group*(1 + 2.0*x[j]);
			g[i + 1] -= 2.0*group*(1 + 2.0*x[i + 1]);
		}
		secGroup -= next;
		secGroup -= x[n - 7] * (1 + x[n - 7]);
		secGroup += current;
		group = x[n - 1] * (2.0 + 5.0* x[n - 1] * x[n - 1]) + 1 - secGroup;
		fx += group * group;
		g[n - 1] += 2.0*group*(2.0 + 15.0 * x[n - 1] * x[n - 1]);

		for (int j = n - 6; j < n - 1; j++)
			g[j] -= 2.0*group*(1 + 2.0*x[j]);

		return fx;
	}
	void InitializeBRYBND
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = -1.0;
	}
	/*========= 4 ========= Chained Wood - CHAINWOO ======== 4000 ==========*/
	double CHAINWOORanged(double *x, double *g, int begin, int end) {
		double fx = 0.0;
		double item1, item2, item3, item4, item5, item6;
		for (int i = begin; i < end; i++)
			g[i] = 0;

		if (begin % 2 != 0) {
			begin++;
		}
		for (int i = begin; i < end - 3; i += 2)
		{
			item1 = (x[i + 1] - x[i] * x[i]);
			item2 = (1 - x[i]);
			item3 = (x[i + 3] - pow(x[i + 2], 2.0));
			item4 = (1 - x[i + 2]);
			item5 = (x[i + 1] + x[i + 3] - 2.0);
			item6 = (x[i + 1] - x[i + 3]);
			fx += 100 * item1*item1 + item2 * item2 + 90 * item3*item3
				+ item4 * item4 + 10.0*item5*item5 + 0.1*item6*item6;
			g[i] -= (400 * item1*x[i] + 2 * item2);
			g[i + 1] += 200 * item1 + 20.0*item5 + 0.2*item6;
			g[i + 2] -= (360 * item3*x[i + 2] + 2.0*item4);
			g[i + 3] += 180 * item3 + 20.0*item5 - 0.2*item6;
		}
		return fx;
	}

	double CHAINWOOHelper(double *x, double *g, int n, int i) {
		if (i % 2 == 1) {
			return 0;
		}
		double fx = 0.0;
		double item1, item2, item3, item4, item5, item6;
		item1 = (x[i + 1] - x[i] * x[i]);
		item2 = (1 - x[i]);
		item3 = (x[i + 3] - pow(x[i + 2], 2.0));
		item4 = (1 - x[i + 2]);
		item5 = (x[i + 1] + x[i + 3] - 2.0);
		item6 = (x[i + 1] - x[i + 3]);
		fx += 100 * item1*item1 + item2 * item2 + 90 * item3*item3
			+ item4 * item4 + 10.0*item5*item5 + 0.1*item6*item6;
		g[i] -= (400 * item1*x[i] + 2 * item2);
		g[i + 1] += 200 * item1 + 20.0*item5 + 0.2*item6;
		g[i + 2] -= (360 * item3*x[i + 2] + 2.0*item4);
		g[i + 3] += 180 * item3 + 20.0*item5 - 0.2*item6;
		return fx;
	}

	double CHAINWOO(double *x, double *g, const int n) {
		double fx = 1.0 + valgradPool(CHAINWOORanged, x, g, n);
		for (int i = block_size - 3; i < n - 3; i = i + block_size) {
			fx += CHAINWOOHelper(x, g, n, i);
			fx += CHAINWOOHelper(x, g, n, i + 1);
			fx += CHAINWOOHelper(x, g, n, i + 2);
		}
		return fx;
	}

	void InitialCHAINWOO
	(
		double *x,
		const int n
	)
	{
		x[0] = -3.0;
		x[1] = -1.0;
		x[2] = -3.0;
		x[3] = -1.0;
		for (int i = 4; i < n; i++)
			x[i] = -2.0;
	}
	/*========= 5 ========== COSINE Function =========== 10000 ===========*/
	double COSINERanged
	(
		double *x,
		double *g,
		int begin,
		int end
	)
	{
		double item;
		double fx(0.0);

		for (int i = begin; i < end; i++)
			g[i] = 0.0;

		for (int i = begin; i < end - 1; i++)
		{
			item = -0.5*x[i + 1] + x[i] * x[i];
			fx += cos(item);
			g[i] -= 2.0*sin(item)*x[i];
			g[i + 1] += 0.5*sin(item);
		}
		return fx;
	}
	double COSINE
	(
		double *x,
		double *g,
		const int n
	)
	{
		double item;
		double fx = valgradPool(COSINERanged, x, g, n);
		for (int i = block_size - 1; i < n - 1; i += block_size)
		{
			item = -0.5*x[i + 1] + x[i] * x[i];
			fx += cos(item);
			g[i] -= 2.0*sin(item)*x[i];
			g[i + 1] += 0.5*sin(item);
		}
		return fx;
	}
	void InitialCOSINE
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = 1.0;
	}
	/*========= 6 ========= CRAGGLVY  =========== 5000 =============*/
	double CRAGGLVYRanged(double *x, double *g, int begin, int end) {
		double fx = 0;
		double item1, item2, element, item3, item4;
		double item1Squared, item2Squared, item3Squared, item4Squared;
		double xipow2, xipow4;
		for (int i = begin; i < end; i++)
			g[i] = 0;

		if (begin % 2 != 0) {
			begin++;
		}
		for (int i = begin; i < end - 3; i += 2)
		{
			item1 = exp(x[i]) - x[i + 1];
			item2 = x[i + 1] - x[i + 2];
			element = x[i + 2] - x[i + 3];
			item3 = tan(element) + element;
			item4 = (x[i + 3] - 1);

			item1Squared = item1 * item1;
			item2Squared = item2 * item2;
			item3Squared = item3 * item3;
			item4Squared = item4 * item4;
			xipow2 = x[i] * x[i];
			xipow4 = xipow2 * xipow2;

			fx += item1Squared * item1Squared + 100 * item2Squared*item2Squared*item2Squared
				+ item3Squared * item3Squared + xipow4 * xipow4 + item4Squared;

			g[i] += 4 * item1Squared*item1*exp(x[i]) + 8 * x[i] * xipow4*xipow2;
			g[i + 1] += -4 * item1*item1Squared + 600 * item2 *item2Squared*item2Squared;
			g[i + 2] += -600 * item2 *item2Squared*item2Squared
				+ 4 * item3* item3Squared*(1 / pow(cos(element), 2.0) + 1);

			g[i + 3] += -4 * item3* item3Squared*(1 / pow(cos(element), 2.0) + 1)
				+ 2 * item4;

		}
		return fx;
	}

	double CRAGGLVYHelper(double *x, double *g, int n, int i) {
		if (i % 2 == 1) {
			return 0;
		}
		double fx = 0.0;
		double item1, item2, element, item3, item4;
		double item1Squared, item2Squared, item3Squared, item4Squared;
		double xipow2, xipow4;

		item1 = exp(x[i]) - x[i + 1];
		item2 = x[i + 1] - x[i + 2];
		element = x[i + 2] - x[i + 3];
		item3 = tan(element) + element;
		item4 = (x[i + 3] - 1);

		item1Squared = item1 * item1;
		item2Squared = item2 * item2;
		item3Squared = item3 * item3;
		item4Squared = item4 * item4;
		xipow2 = x[i] * x[i];
		xipow4 = xipow2 * xipow2;

		fx += item1Squared * item1Squared + 100 * item2Squared*item2Squared*item2Squared
			+ item3Squared * item3Squared + xipow4 * xipow4 + item4Squared;

		g[i] += 4 * item1Squared*item1*exp(x[i]) + 8 * x[i] * xipow4*xipow2;
		g[i + 1] += -4 * item1*item1Squared + 600 * item2 *item2Squared*item2Squared;
		g[i + 2] += -600 * item2 *item2Squared*item2Squared
			+ 4 * item3* item3Squared*(1 / pow(cos(element), 2.0) + 1);

		g[i + 3] += -4 * item3* item3Squared*(1 / pow(cos(element), 2.0) + 1)
			+ 2 * item4;
		return fx;
	}

	double CRAGGLVY(double *x, double *g, const int n) {
		double fx = 1.0 + valgradPool(CRAGGLVYRanged, x, g, n);
		for (int i = block_size - 3; i < n - 3; i = i + block_size) {
			fx += CRAGGLVYHelper(x, g, n, i);
			fx += CRAGGLVYHelper(x, g, n, i + 1);
			fx += CRAGGLVYHelper(x, g, n, i + 2);
		}
		return fx;
	}

	void InitialCRAGGLVY
	(
		double *x,
		const int n
	)
	{
		x[0] = 1.0;

		for (int i = 1; i < n; i++)
			x[i] = 2.0;
	}


	/*===== 7 =============== CURLY10 Function ============= 500 ===================*/
	double CURLY10
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = 0.0;
		double cube;
		int k(10);
		int i, j;
		double q;
		for (int i = 0; i < n; i++)
			g[i] = 0;

		q = 0.0;
		for (int j = 0; j <= k; j++)
			q += x[j];
		fx += pow(q, 4.0) - 20 * q*q - 0.1*q;
		for (int j = 0; j <= k; j++)
			g[j] += 4 * pow(q, 3.0) - 40 * q - 0.1;

		for (i = 1; i < n - k; i++)
		{
			q = q - x[i - 1] + x[i + k];
			double t0 = q, t1 = t0 * t0, t2 = t1 * t1;
			fx += t2 - 20 * t1 - 0.1*t0;;
			cube = 4 * t0*t1 - 40 * t0 - 0.1;
			j = i;
			g[j] += cube; ++j;
			for (; j <= i + k; )
			{
				g[j] += cube; ++j;
				g[j] += cube; ++j;
				g[j] += cube; ++j;
				g[j] += cube; ++j;
				g[j] += cube; ++j;
			}
		}
		for (; i < n; i++)
		{
			q -= x[i - 1];
			double t0 = q, t1 = t0 * t0, t2 = t1 * t1;
			fx += t2 - 20 * t1 - 0.1*t0;;
			cube = 4 * t0*t1 - 40 * t0 - 0.1;
			for (int j = i; j <= n - 1; j++)
				g[j] += cube;
		}
		return fx;
	}
	void InitializeCURLY
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = 0.0001 / (n + 1);
	}
	/*===== 8 ============= CURLY20 Function ============= 500 ===================*/
	double CURLY20
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = 0.0;
		double cube;
		int k(20);
		int i, j;
		double q;
		for (int i = 0; i < n; i++)
			g[i] = 0;

		q = 0.0;
		for (int j = 0; j <= k; j++)
			q += x[j];
		fx += pow(q, 4.0) - 20 * q*q - 0.1*q;
		for (int j = 0; j <= k; j++)
			g[j] += 4 * pow(q, 3.0) - 40 * q - 0.1;

		for (i = 1; i < n - k; i++)
		{
			q = q - x[i - 1] + x[i + k];
			double t0 = q, t1 = t0 * t0, t2 = t1 * t1;
			fx += t2 - 20 * t1 - 0.1*t0;;
			cube = 4 * t0*t1 - 40 * t0 - 0.1;
			j = i;
			g[j] += cube; ++j;
			for (; j <= i + k; )
			{
				g[j] += cube; ++j;
				g[j] += cube; ++j;
				g[j] += cube; ++j;
				g[j] += cube; ++j;
				g[j] += cube; ++j;
			}
		}
		for (; i < n; i++)
		{
			q -= x[i - 1];
			double t0 = q, t1 = t0 * t0, t2 = t1 * t1;
			fx += t2 - 20 * t1 - 0.1*t0;;
			cube = 4 * t0*t1 - 40 * t0 - 0.1;
			for (int j = i; j <= n - 1; j++)
				g[j] += cube;
		}
		return fx;
	}
	//initialization in curly functions is common
	/*===== 9 ============= CURLY30 Function ======== 500 =================*/
	double CURLY30
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = 0.0;
		double cube;
		int k(30);
		int i, j;
		double q;
		for (int i = 0; i < n; i++)
			g[i] = 0;

		q = 0.0;
		for (int j = 0; j <= k; j++)
			q += x[j];
		fx += pow(q, 4.0) - 20 * q*q - 0.1*q;
		cube = 4 * pow(q, 3.0) - 40 * q - 0.1;
		for (j = 0; j <= k; j++)
			g[j] += cube;

		for (i = 1; i < n - k; i++)
		{
			q = q - x[i - 1] + x[i + k];
			double t0 = q, t1 = t0 * t0, t2 = t1 * t1;
			fx += t2 - 20 * t1 - 0.1*t0;;
			cube = 4 * t0*t1 - 40 * t0 - 0.1;
			j = i;
			g[j] += cube; ++j;
			for (; j <= i + k; )
			{
				g[j] += cube; ++j;
				g[j] += cube; ++j;
				g[j] += cube; ++j;
				g[j] += cube; ++j;
				g[j] += cube; ++j;
			}
		}
		for (; i < n; i++)
		{
			q -= x[i - 1];
			double t0 = q, t1 = t0 * t0, t2 = t1 * t1;
			fx += t2 - 20 * t1 - 0.1*t0;;
			cube = 4 * t0*t1 - 40 * t0 - 0.1;
			for (j = i; j <= n - 1; j++)
				g[j] += cube;
		}
		return fx;
	}

	/*===== 10 ======================= DIXMAANA Function ================== 5000 ================*/
	double DIXMAANARanged1
	(
		double *x,
		double *g,
		int begin,
		int end
	)
	{
		for (int i = begin; i < end; i++)
			g[i] = 0;
		return 0;
	}
	double DIXMAANARanged2
	(
		double *x,
		double *g,
		int begin,
		int end,
		int m,
		int placeholder
	)
	{
		double fx = 0.0;
		double item(0.0);
		for (int i = begin; i < end; i++)
		{
			item = x[i] * pow(x[i + m], 2);
			fx += pow(x[i], 2) + 0.125*pow(item, 2) + 0.125*x[i] * x[i + 2 * m];
			g[i] += 2 * x[i] + 0.25*item *pow(x[i + m], 2) + 0.125*x[i + 2 * m];
			g[i + m] += 0.5*item *x[i] * x[i + m];
			g[i + 2 * m] += 0.125*x[i];
		}
		return fx;
	}
	double DIXMAANARanged3
	(
		double *x,
		double *g,
		int begin,
		int end,
		int m,
		int placeholder
	)
	{
		double fx = 0.0;
		double item(0.0);
		for (int i = begin; i < end; i++)
		{
			item = x[i] * pow(x[i + m], 2);
			fx += pow(x[i], 2) + 0.125*pow(item, 2);
			g[i] += 2 * x[i] + 0.25*item *pow(x[i + m], 2);
			g[i + m] += 0.5*item *x[i] * x[i + m];
		}
		return fx;
	}
	double DIXMAANARanged4
	(
		double *x,
		double *g,
		int begin,
		int end,
		int m,
		int placeholder
	)
	{
		double fx = 0.0;
		double item(0.0);
		for (int i = begin; i < end; i++)
		{
			fx += pow(x[i], 2);
			g[i] += 2 * x[i];
		}
		return fx;
	}
	double DIXMAANA
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = 1.0;
		double item(0.0);
		int m(n / 3);
		initializer(n, 100);
		valgradPool(DIXMAANARanged1, x, g, n);
		initializer(m, 100);
		fx += valgradPoolRangedWithM(DIXMAANARanged2, x, g, 0, m, m);
		initializer(m, 100);
		fx += valgradPoolRangedWithM(DIXMAANARanged3, x, g, m, 2*m, m);
		initializer(n - (2*m), 100);
		fx += valgradPoolRangedWithM(DIXMAANARanged4, x, g, 2 * m, n, m);
		return fx;
	}
	void InitializeDIXMAAN
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = 2.0;
	}

	/*===== 11 ======================= DIXMAANB Function ================== 5000 ================*/
	double DIXMAANBRanged1
	(
		double *x,
		double *g,
		int begin,
		int end
	)
	{
		for (int i = begin; i < end; i++)
			g[i] = 0;
		return 0;
	}
	double DIXMAANBHelper2(double *x, double *g, int i, int m)
	{
		double item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		double item2 = 0.25* x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.5*item2*pow(x[i + m], 2) +
			(0.0625*x[i + 2 * m]);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 * x[i] * x[i + m];
		g[i + 2 * m] += 0.0625*x[i];
		return pow(x[i], 2) + item1 * item1 + item2 * item2 + (0.0625*x[i] * x[i + 2 * m]);
	}
	double DIXMAANBRanged2
	(
		double *x,
		double *g,
		int begin,
		int end,
		int m,
		int placeholder
	)
	{
		double fx = 0.0;
		for (int i = begin; i < end - 1; i++)
		{
			fx += DIXMAANBHelper2(x, g, i, m);
		}
		return fx;
	}
	double DIXMAANBHelper3(double *x, double *g, int i, int m)
	{
		double item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		double item2 = 0.25* x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.5*item2*pow(x[i + m], 2);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 * x[i] * x[i + m];
		return pow(x[i], 2) + item1 * item1 + item2 * item2;
	}
	double DIXMAANBRanged3
	(
		double *x,
		double *g,
		int begin,
		int end,
		int m,
		int placeholder
	)
	{
		double fx = 0.0;
		for (int i = begin; i < end - 1; i++)
		{
			fx += DIXMAANBHelper3(x, g, i, m);
		}
		return fx;
	}
	double DIXMAANBHelper4(double *x, double *g, int i, int m)
	{
		double item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i] += 2 * x[i] + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		return pow(x[i], 2) + item1 * item1;
	}
	double DIXMAANBRanged4
	(
		double *x,
		double *g,
		int begin,
		int end,
		int m,
		int placeholder
	)
	{
		double fx = 0.0;
		for (int i = begin; i < end - 1; i++)
		{
			fx += DIXMAANBHelper4(x, g, i, m);
		}
		return fx;
	}
	double DIXMAANB
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = 1.0;
		double item1(0.0);
		double item2(0.0);
		int m(n / 3);
		int skippedIndex;

		initializer(n, 100);
		valgradPool(DIXMAANBRanged1, x, g, n);
		initializer(m, 100);
		fx += valgradPoolRangedWithM(DIXMAANBRanged2, x, g, 0, m, m);
		initializer(m, 100);
		fx += valgradPoolRangedWithM(DIXMAANBRanged3, x, g, m, 2 * m, m);
		initializer(n - (2 * m), 100);
		fx += valgradPoolRangedWithM(DIXMAANBRanged4, x, g, 2 * m, n, m);

		initializer(m, 100);
		skippedIndex = block_size - 1;
		for (int i = 0; i < (num_threads - 1); i++)
		{
			fx += DIXMAANBHelper2(x, g, skippedIndex, m);
			skippedIndex += block_size;
		}
		fx += DIXMAANBHelper2(x, g, m - 1, m);

		initializer(m, 100);
		skippedIndex = m + block_size - 1;
		for (int i = 0; i < (num_threads - 1); i++)
		{
			fx += DIXMAANBHelper3(x, g, skippedIndex, m);
			skippedIndex += block_size;
		}
		fx += DIXMAANBHelper3(x, g, (2 * m) - 1, m);

		initializer(n - (2 * m), 100);
		skippedIndex = (2 * m) + block_size - 1;
		for (int i = 0; i < (num_threads - 1); i++)
		{
			fx += DIXMAANBHelper4(x, g, skippedIndex, m);
			skippedIndex += block_size;
		}

		fx += pow(x[n - 1], 2);
		g[n - 1] += 2 * x[n - 1];
		return fx;
	}
	/*=====12  ======================= DIXMAANC Function ================== 5000 ================*/
	double DIXMAANC
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = 1.0;
		double item1(0.0);
		double item2(0.0);
		int m(n / 3);
		for (int i = 0; i < n; i++)
			g[i] = 0;
		for (int i = 0; i < m; i++)
		{
			item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			item2 = x[i] * pow(x[i + m], 2);
			fx += pow(x[i], 2) + 0.125*item1*item1 + 0.125*item2*item2 + 0.125*x[i] * x[i + 2 * m];
			g[i] += 2 * x[i] + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.25*item2*pow(x[i + m], 2) + 0.125*x[i + 2 * m];
			g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
			g[i + m] += 0.5* item2 *x[i] * x[i + m];
			g[i + 2 * m] += 0.125*x[i];
		}
		for (int i = m; i < 2 * m; i++)
		{
			item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			item2 = x[i] * pow(x[i + m], 2);
			fx += pow(x[i], 2) + 0.125*item1*item1 + 0.125*item2*item2;
			g[i] += 2 * x[i] + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.25*item2*pow(x[i + m], 2);
			g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
			g[i + m] += 0.5* item2 *x[i] * x[i + m];
		}
		for (int i = 2 * m; i < n - 1; i++)
		{
			item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			fx += pow(x[i], 2) + 0.125*item1*item1;
			g[i] += 2 * x[i] + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
			g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
		}
		fx += pow(x[n - 1], 2);
		g[n - 1] += 2 * x[n - 1];
		return fx;
	}
	/*===== 13 ======================= DIXMAAND Function ================== 5000 ================*/
	double DIXMAAND
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = 1.0;
		double item1(0.0);
		double item2(0.0);
		int m(n / 3);
		for (int i = 0; i < n; i++)
			g[i] = 0;
		for (int i = 0; i < m; i++)
		{
			item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			item2 = x[i] * pow(x[i + m], 2);
			fx += pow(x[i], 2) + 0.26*item1*item1 + 0.26*item2*item2 + 0.26*x[i] * x[i + 2 * m];
			g[i] += 2 * x[i] + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.52*item2*pow(x[i + m], 2) + 0.26*x[i + 2 * m];
			g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
			g[i + m] += 1.04* item2 *x[i] * x[i + m];
			g[i + 2 * m] += 0.26*x[i];
		}
		for (int i = m; i < 2 * m; i++)
		{
			item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			item2 = x[i] * pow(x[i + m], 2);
			fx += pow(x[i], 2) + 0.26*item1*item1 + 0.26*item2*item2;
			g[i] += 2 * x[i] + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.52*item2*pow(x[i + m], 2);
			g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
			g[i + m] += 1.04* item2 *x[i] * x[i + m];
		}
		for (int i = 2 * m; i < n - 1; i++)
		{
			item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			fx += pow(x[i], 2) + 0.26*item1*item1;
			g[i] += 2 * x[i] + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
			g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
		}
		fx += pow(x[n - 1], 2);
		g[n - 1] += 2 * x[n - 1];
		return fx;
	}
	/*===== 14 ======================= DIXMAANE Function ================== 3000 ================*/
	double DIXMAANE
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = 1.0;
		double item(0.0);
		int m(n / 3);
		for (int i = 0; i < n; i++)
			g[i] = 0;
		for (int i = 0; i < m; i++)
		{
			item = x[i] * pow(x[i + m], 2);
			fx += (pow(x[i], 2)*(i + 1)) / n + 0.125*pow(item, 2) + (0.125*x[i] * x[i + 2 * m] * (i + 1)) / n;
			g[i] += (2 * x[i] * (i + 1)) / n + 0.25*item *pow(x[i + m], 2) + (0.125*x[i + 2 * m] * (i + 1)) / n;
			g[i + m] += 0.5*item *x[i] * x[i + m];
			g[i + 2 * m] += (0.125*x[i] * (i + 1)) / n;
		}
		for (int i = m; i < 2 * m; i++)
		{
			item = x[i] * pow(x[i + m], 2);
			fx += (pow(x[i], 2)*(i + 1)) / n + 0.125*pow(item, 2);
			g[i] += (2 * x[i] * (i + 1)) / n + 0.25*item *pow(x[i + m], 2);
			g[i + m] += 0.5*item *x[i] * x[i + m];
		}
		for (int i = 2 * m; i < n; i++)
		{
			fx += (pow(x[i], 2)*i) / n;
			g[i] += (2 * x[i] * (i + 1)) / n;
		}
		return fx;
	}
	/*===== 15 ======================= DIXMAANF Function ================== 3000 ================*/
	double DIXMAANF
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = 1.0;
		double item1(0.0);
		double item2(0.0);
		int m(n / 3);
		for (int i = 0; i < n; i++)
			g[i] = 0;
		for (int i = 0; i < m; i++)
		{
			item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			item2 = 0.25* x[i] * pow(x[i + m], 2);
			fx += (pow(x[i], 2)*(i + 1)) / n + item1 * item1 + item2 * item2 + (0.0625*x[i] * x[i + 2 * m] * (i + 1)) / n;
			g[i] += (2 * x[i] * (i + 1)) / n + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.5*item2*pow(x[i + m], 2) +
				(0.0625*x[i + 2 * m] * (i + 1)) / n;
			g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
			g[i + m] += item2 * x[i] * x[i + m];
			g[i + 2 * m] += (0.0625*x[i] * (i + 1)) / n;
		}
		for (int i = m; i < 2 * m; i++)
		{
			item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			item2 = 0.25* x[i] * pow(x[i + m], 2);
			fx += (pow(x[i], 2)*(i + 1)) / n + item1 * item1 + item2 * item2;
			g[i] += (2 * x[i] * (i + 1)) / n + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.5*item2*pow(x[i + m], 2);
			g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
			g[i + m] += item2 * x[i] * x[i + m];
		}
		for (int i = 2 * m; i < n - 1; i++)
		{
			item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			fx += (pow(x[i], 2)*(i + 1)) / n + item1 * item1;
			g[i] += (2 * x[i] * (i + 1)) / n + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
			g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		}
		fx += pow(x[n - 1], 2);
		g[n - 1] += 2 * x[n - 1];
		return fx;
	}
	/*===== 16 ======================= DIXMAANG Function ================== 5000 ================*/
	double DIXMAANG
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = 1.0;
		double item1(0.0);
		double item2(0.0);
		int m(n / 3);
		for (int i = 0; i < n; i++)
			g[i] = 0;
		for (int i = 0; i < m; i++)
		{
			item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			item2 = x[i] * pow(x[i + m], 2);
			fx += (pow(x[i], 2)*(i + 1)) / n + 0.125*item1*item1 + 0.125*item2*item2
				+ (0.125*x[i] * x[i + 2 * m] * (i + 1)) / n;
			g[i] += (2 * x[i] * (i + 1)) / n + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.25*item2*pow(x[i + m], 2)
				+ (0.125*x[i + 2 * m] * (i + 1)) / n;
			g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
			g[i + m] += 0.5* item2 *x[i] * x[i + m];
			g[i + 2 * m] += (0.125*x[i] * (i + 1)) / n;
		}
		for (int i = m; i < 2 * m; i++)
		{
			item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			item2 = x[i] * pow(x[i + m], 2);
			fx += (pow(x[i], 2)*(i + 1)) / n + 0.125*item1*item1 + 0.125*item2*item2;
			g[i] += (2 * x[i] * (i + 1)) / n + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.25*item2*pow(x[i + m], 2);
			g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
			g[i + m] += 0.5* item2 *x[i] * x[i + m];
		}
		for (int i = 2 * m; i < n - 1; i++)
		{
			item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			fx += (pow(x[i], 2)*(i + 1)) / n + 0.125*item1*item1;
			g[i] += (2 * x[i] * (i + 1)) / n + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
			g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
		}
		fx += pow(x[n - 1], 2);
		g[n - 1] += 2 * x[n - 1];
		return fx;
	}

	/*===== 17 ======================= DIXMAANH Function ================== 5000 ================*/
	double DIXMAANH
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = 1.0;
		double item1(0.0);
		double item2(0.0);
		int m(n / 3);
		for (int i = 0; i < n; i++)
			g[i] = 0;
		for (int i = 0; i < m; i++)
		{
			item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			item2 = x[i] * pow(x[i + m], 2);
			fx += (pow(x[i], 2)*(i + 1)) / n + 0.26*item1*item1 + 0.26*item2*item2
				+ (0.26*x[i] * x[i + 2 * m] * (i + 1)) / n;
			g[i] += (2 * x[i] * (i + 1)) / n + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.52*item2*pow(x[i + m], 2)
				+ (0.26*x[i + 2 * m] * (i + 1)) / n;
			g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
			g[i + m] += 1.04* item2 *x[i] * x[i + m];
			g[i + 2 * m] += (0.26*x[i] * (i + 1)) / n;
		}
		for (int i = m; i < 2 * m; i++)
		{
			item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			item2 = x[i] * pow(x[i + m], 2);
			fx += (pow(x[i], 2)*(i + 1)) / n + 0.26*item1*item1 + 0.26*item2*item2;
			g[i] += (2 * x[i] * (i + 1)) / n + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]) + 0.52*item2*pow(x[i + m], 2);
			g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
			g[i + m] += 1.04* item2 *x[i] * x[i + m];
		}
		for (int i = 2 * m; i < n - 1; i++)
		{
			item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			fx += (pow(x[i], 2)*(i + 1)) / n + 0.26*item1*item1;
			g[i] += (2 * x[i] * (i + 1)) / n + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
			g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
		}
		fx += pow(x[n - 1], 2);
		g[n - 1] += 2 * x[n - 1];
		return fx;
	}
	/*===== 18 ======================= DIXMAANI Function ================== 3000 ================*/
	double DIXMAANIRanged1
	(
		double *x,
		double *g,
		int begin,
		int end
	)
	{
		for (int i = begin; i < end; i++)
			g[i] = 0;
		return 0;
	}
	double DIXMAANIRanged2
	(
		double *x,
		double *g,
		int begin,
		int end,
		int m,
		int n
	)
	{
		double fx = 0.0;
		double item(0.0);
		for (int i = begin; i < end; i++)
		{
			item = x[i] * pow(x[i + m], 2);
			fx += pow(x[i], 2)*pow((i + 1) / double(n), 2) + 0.125*pow(item, 2)
				+ 0.125*x[i] * x[i + 2 * m] * pow((i + 1) / double(n), 2);
			g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.25*item *pow(x[i + m], 2)
				+ 0.125*x[i + 2 * m] * pow((i + 1) / double(n), 2);
			g[i + m] += 0.5*item *x[i] * x[i + m];
			g[i + 2 * m] += 0.125*x[i] * pow((i + 1) / double(n), 2);
		}
		return fx;
	}
	double DIXMAANIRanged3
	(
		double *x,
		double *g,
		int begin,
		int end,
		int m,
		int n
	)
	{
		double fx = 0.0;
		double item(0.0);
		for (int i = begin; i < end; i++)
		{
			item = x[i] * pow(x[i + m], 2);
			fx += pow(x[i], 2)*pow((i + 1) / double(n), 2) + 0.125*pow(item, 2);
			g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.25*item *pow(x[i + m], 2);
			g[i + m] += 0.5*item *x[i] * x[i + m];
		}
		return fx;
	}
	double DIXMAANIRanged4
	(
		double *x,
		double *g,
		int begin,
		int end,
		int m,
		int n
	)
	{
		double fx = 0.0;
		double item(0.0);
		for (int i = begin; i < end; i++)
		{
			fx += pow(x[i], 2)*pow((i + 1) / double(n), 2);
			g[i] += 2 * x[i] * pow((i + 1) / double(n), 2);
		}
		return fx;
	}

	double DIXMAANI
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = 1.0;
		double item(0.0);
		int m(n / 3);

		initializer(n, 100);
		valgradPool(DIXMAANIRanged1, x, g, n);
		initializer(m, 100);
		fx += valgradPoolRangedWithMWithN(DIXMAANIRanged2, x, g, 0, m, m, n);
		initializer(m, 100);
		fx += valgradPoolRangedWithMWithN(DIXMAANIRanged3, x, g, m, 2 * m, m, n);
		initializer(n - (2 * m), 100);
		fx += valgradPoolRangedWithMWithN(DIXMAANIRanged4, x, g, 2 * m, n, m, n);
		return fx;
	}
	/*===== 19 ======================= DIXMAANJ Function ================== 3000 ================*/
	double DIXMAANJRanged1
	(
		double *x,
		double *g,
		int begin,
		int end
	)
	{
		for (int i = begin; i < end; i++)
			g[i] = 0;
		return 0;
	}
	double DIXMAANJHelper2(double *x, double *g, int i, int m, int n)
	{
		double item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		double item2 = 0.25* x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.5*item2*pow(x[i + m], 2) + 0.0625*x[i + 2 * m] * pow((i + 1) / double(n), 2);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 * x[i] * x[i + m];
		g[i + 2 * m] += 0.0625*x[i] * pow((i + 1) / double(n), 2);
		return pow(x[i], 2)*pow((i + 1) / double(n), 2) + item1 * item1 + item2 * item2
			+ 0.0625*x[i] * x[i + 2 * m] * pow((i + 1) / double(n), 2);
	}
	double DIXMAANJRanged2
	(
		double *x,
		double *g,
		int begin,
		int end,
		int m,
		int n
	)
	{
		double fx = 0.0;
		for (int i = begin; i < end - 1; i++)
		{
			fx += DIXMAANJHelper2(x, g, i, m, n);
		}
		return fx;
	}
	double DIXMAANJHelper3(double *x, double *g, int i, int m, int n)
	{
		double item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		double item2 = 0.25* x[i] * pow(x[i + m], 2);
		g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1])
			+ 0.5*item2*pow(x[i + m], 2);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		g[i + m] += item2 * x[i] * x[i + m];
		return pow(x[i], 2)*pow((i + 1) / double(n), 2) + item1 * item1 + item2 * item2;
	}
	double DIXMAANJRanged3
	(
		double *x,
		double *g,
		int begin,
		int end,
		int m,
		int n
	)
	{
		double fx = 0.0;
		for (int i = begin; i < end - 1; i++)
		{
			fx += DIXMAANJHelper3(x, g, i, m, n);
		}
		return fx;
	}
	double DIXMAANJHelper4(double *x, double *g, int i, int m, int n)
	{
		double item1 = 0.25* x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
		g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.5*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
		g[i + 1] += 0.5*item1*x[i] * (1 + 2 * x[i + 1]);
		return pow(x[i], 2)*pow((i + 1) / double(n), 2) + item1 * item1;
	}
	double DIXMAANJRanged4
	(
		double *x,
		double *g,
		int begin,
		int end,
		int m,
		int n
	)
	{
		double fx = 0.0;
		for (int i = begin; i < end - 1; i++)
		{
			fx += DIXMAANJHelper4(x, g, i, m, n);
		}
		return fx;
	}

	double DIXMAANJ
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = 1.0;
		double item1(0.0);
		double item2(0.0);
		int m(n / 3);
		int skippedIndex;

		initializer(n, 100);
		valgradPool(DIXMAANJRanged1, x, g, n);
		initializer(m, 100);
		fx += valgradPoolRangedWithMWithN(DIXMAANJRanged2, x, g, 0, m, m, n);
		initializer(m, 100);
		fx += valgradPoolRangedWithMWithN(DIXMAANJRanged3, x, g, m, 2 * m, m, n);
		initializer(n - (2 * m), 100);
		fx += valgradPoolRangedWithMWithN(DIXMAANJRanged4, x, g, 2 * m, n, m, n);

		initializer(m, 100);
		skippedIndex = block_size - 1;
		for (int i = 0; i < (num_threads - 1); i++)
		{
			fx += DIXMAANJHelper2(x, g, skippedIndex, m, n);
			skippedIndex += block_size;
		}
		fx += DIXMAANJHelper2(x, g, m - 1, m, n);

		initializer(m, 100);
		skippedIndex = m + block_size - 1;
		for (int i = 0; i < (num_threads - 1); i++)
		{
			fx += DIXMAANJHelper3(x, g, skippedIndex, m, n);
			skippedIndex += block_size;
		}
		fx += DIXMAANJHelper3(x, g, (2 * m) - 1, m, n);

		initializer(n - (2 * m), 100);
		skippedIndex = (2 * m) + block_size - 1;
		for (int i = 0; i < (num_threads - 1); i++)
		{
			fx += DIXMAANJHelper4(x, g, skippedIndex, m, n);
			skippedIndex += block_size;
		}


		fx += pow(x[n - 1], 2);
		g[n - 1] += 2 * x[n - 1];
		return fx;
	}
	/*===== 20 ======================= DIXMAANK Function ================== 5000 ================*/
	double DIXMAANK
	(

		double *x,
		double *g,
		const int n
	)
	{
		double fx = 1.0;
		double item1(0.0);
		double item2(0.0);
		int m(n / 3);
		for (int i = 0; i < n; i++)
			g[i] = 0;
		for (int i = 0; i < m; i++)
		{
			item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			item2 = x[i] * pow(x[i + m], 2);
			fx += pow(x[i], 2)*pow((i + 1) / double(n), 2) + 0.125*item1*item1 + 0.125*item2*item2
				+ 0.125*x[i] * x[i + 2 * m] * pow((i + 1) / double(n), 2);
			g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1])
				+ 0.25*item2*pow(x[i + m], 2) + 0.125*x[i + 2 * m] * pow((i + 1) / double(n), 2);
			g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
			g[i + m] += 0.5* item2 *x[i] * x[i + m];
			g[i + 2 * m] += 0.125*x[i] * pow((i + 1) / double(n), 2);
		}
		for (int i = m; i < 2 * m; i++)
		{
			item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			item2 = x[i] * pow(x[i + m], 2);
			fx += pow(x[i], 2)*pow((i + 1) / double(n), 2) + 0.125*item1*item1 + 0.125*item2*item2;
			g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1])
				+ 0.25*item2*pow(x[i + m], 2);
			g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
			g[i + m] += 0.5* item2 *x[i] * x[i + m];
		}
		for (int i = 2 * m; i < n - 1; i++)
		{
			item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			fx += pow(x[i], 2)*pow((i + 1) / double(n), 2) + 0.125*item1*item1;
			g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.25*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
			g[i + 1] += 0.25*item1*x[i] * (1 + 2 * x[i + 1]);
		}
		fx += pow(x[n - 1], 2);
		g[n - 1] += 2 * x[n - 1];
		return fx;
	}

	/*===== 21 =================== DIXMAANL Function ================== 5000 ================*/
	double DIXMAANL
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = 1.0;
		double item1(0.0);
		double item2(0.0);
		int m(n / 3);
		for (int i = 0; i < n; i++)
			g[i] = 0;
		for (int i = 0; i < m; i++)
		{
			item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			item2 = x[i] * pow(x[i + m], 2);
			fx += pow(x[i], 2)*pow((i + 1) / double(n), 2) + 0.26*item1*item1 + 0.26*item2*item2
				+ 0.26*x[i] * x[i + 2 * m] * pow((i + 1) / double(n), 2);
			g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1])
				+ 0.52*item2*pow(x[i + m], 2) + 0.26*x[i + 2 * m] * pow((i + 1) / double(n), 2);
			g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
			g[i + m] += 1.04* item2 *x[i] * x[i + m];
			g[i + 2 * m] += 0.26*x[i] * pow((i + 1) / double(n), 2);
		}
		for (int i = m; i < 2 * m; i++)
		{
			item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			item2 = x[i] * pow(x[i + m], 2);
			fx += pow(x[i], 2)*pow((i + 1) / double(n), 2) + 0.26*item1*item1 + 0.26*item2*item2;
			g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1])
				+ 0.52*item2*pow(x[i + m], 2);
			g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
			g[i + m] += 1.04* item2 *x[i] * x[i + m];
		}
		for (int i = 2 * m; i < n - 1; i++)
		{
			item1 = x[i] * (x[i + 1] + x[i + 1] * x[i + 1]);
			fx += pow(x[i], 2)*pow((i + 1) / double(n), 2) + 0.26*item1*item1;
			g[i] += 2 * x[i] * pow((i + 1) / double(n), 2) + 0.52*item1*(x[i + 1] + x[i + 1] * x[i + 1]);
			g[i + 1] += 0.52*item1*x[i] * (1 + 2 * x[i + 1]);
		}
		fx += pow(x[n - 1], 2);
		g[n - 1] += 2 * x[n - 1];
		return fx;
	}
	/*========= 22 ========== DIXON3DQ Function ============= 5000 ==============*/
	double DIXON3DQRanged
	(
		double *x,
		double *g,
		int begin,
		int end
	)
	{
		double item;
		double fx(0.0);


		for (int i = begin; i < end; i++)
			g[i] = 0.0;

		for (int i = begin + 1; i < end - 1; i++)
		{
			item = (x[i] - x[i + 1]);
			fx += item * item;
			g[i] += 2.0*item;
			g[i + 1] -= 2.0*item;
		}
		return fx;
	}
	double DIXON3DQ
	(
		double *x,
		double *g,
		const int n
	)
	{
		double item;
		double fx = valgradPool(DIXON3DQRanged, x, g, n);

		item = x[0] - 1;
		fx += item * item;
		g[0] += 2.0*item;

		item = x[n - 1] - 1;
		fx += item * item;
		g[n - 1] += 2.0*item;

		for (int i = block_size - 1; i < n - 1; i += block_size)
		{
			item = (x[i] - x[i + 1]);
			fx += item * item;
			g[i] += 2.0*item;
			g[i + 1] -= 2.0*item;
		}
		for (int i = block_size; i < n - 1; i += block_size)
		{
			item = (x[i] - x[i + 1]);
			fx += item * item;
			g[i] += 2.0*item;
			g[i + 1] -= 2.0*item;
		}
		return fx;
	}
	void InitialDIXON3DQ
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = -1.0;
	}
	/*===== 23 ===================== DQDRTIC Function ===== 10000 ======================*/
	double DQDRTICRanged
	(
		double *x,
		double *g,
		int begin,
		int end
	)
	{
		for (int i = begin; i < end; i++)
			g[i] = 0;
		double t0 = x[begin] * x[begin], t1 = x[begin + 1] * x[begin + 1], t2 = x[begin + 2] * x[begin + 2];
		double fx = t0 + 100 * (t1 + t2);
		g[begin] += 2 * x[begin];
		g[begin + 1] += 200 * x[begin + 1];
		g[begin + 2] += 200 * x[begin + 2];
		for (int i = begin + 1; i < end - 2; i++)
		{
			t0 = t1;
			t1 = t2;
			t2 = x[i + 2] * x[i + 2];
			fx += t0 + 100 * (t1 + t2);
			g[i] += 2 * x[i];
			g[i + 1] += 200 * x[i + 1];
			g[i + 2] += 200 * x[i + 2];
		}
		return fx;
	}
	double DQDRTIC
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = valgradPool(DQDRTICRanged, x, g, n);
		double t0, t1, t2;
		for (int i = block_size - 2; i < n - 2; i += block_size)
		{
			t0 = x[i] * x[i];
			t1 = x[i + 1] * x[i + 1];
			t2 = x[i + 2] * x[i + 2];
			fx += t0 + 100 * (t1 + t2);
			g[i] += 2 * x[i];
			g[i + 1] += 200 * x[i + 1];
			g[i + 2] += 200 * x[i + 2];
		}
		for (int i = block_size - 1; i < n - 2; i += block_size)
		{
			t0 = x[i] * x[i];
			t1 = x[i + 1] * x[i + 1];
			t2 = x[i + 2] * x[i + 2];
			fx += t0 + 100 * (t1 + t2);
			g[i] += 2 * x[i];
			g[i + 1] += 200 * x[i + 1];
			g[i + 2] += 200 * x[i + 2];
		}
		return fx;
	}
	void InitialDQDRTIC
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = 3.0;
	}

	/*===== 24 ===================== DQRTIC Function ===== 5000 ======================*/
	double DQRTICRanged(double *x, double *g, int begin, int end) {
		double fx = 0.0;
		double item, squaredItem;
		for (int i = begin; i < end; i++)
		{
			item = x[i] - i - 1;
			squaredItem = item * item;
			fx += squaredItem * squaredItem;
			g[i] = 4 * squaredItem *item;
		}
		return fx;
	}
	double DQRTIC(double *x, double *g, const int n) {
		return valgradPool(DQRTICRanged, x, g, n);
	}
	void InitialDQRTIC
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = 2.0;
	}
	/*===== 25 ===================== EDENSCH Function ===== 2000 ======================*/
	double EDENSCHRanged
	(
		double *x,
		double *g,
		int begin,
		int end
	)
	{
		double fx = 0.0;
		for (int i = begin; i < end; i++)
			g[i] = 0;
		double item1, item2, item3;
		double squaredItem1;
		for (int i = begin; i < end - 1; i++)
		{
			item1 = x[i] - 2;
			item2 = x[i] * x[i + 1] - 2 * x[i + 1];
			item3 = x[i + 1] + 1;
			squaredItem1 = item1 * item1;
			fx += 16 + squaredItem1 * squaredItem1 + item2 * item2 + item3 * item3;
			g[i] += 4 * squaredItem1 * item1 + 2 * item2*x[i + 1];
			g[i + 1] += 2 * item2*(x[i] - 2.0) + 2 * item3;
		}
		return fx;
	}
	double EDENSCH
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = valgradPool(EDENSCHRanged, x, g, n);
		double item1, item2, item3;
		double squaredItem1;
		for (int i = block_size - 1; i < n - 1; i += block_size)
		{
			item1 = x[i] - 2;
			item2 = x[i] * x[i + 1] - 2 * x[i + 1];
			item3 = x[i + 1] + 1;
			squaredItem1 = item1 * item1;
			fx += 16 + squaredItem1 * squaredItem1 + item2 * item2 + item3 * item3;
			g[i] += 4 * squaredItem1 * item1 + 2 * item2*x[i + 1];
			g[i + 1] += 2 * item2*(x[i] - 2.0) + 2 * item3;
		}
		return fx;
	}
	void InitialEDENSCH
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = 0.0;
	}
	/*===== 26 ===================== EG2 Function ===== 5000 ======================*/
	double EG2Ranged(double *x, double *g, double *sum, INT index, INT begin, INT end)
	{
		for (int i = begin; i < end; i++)
			g[i] = 0;
		double fx = 0;
		double item;

		for (int i = begin; i < end - 1; i++)
		{
			item = x[index] + x[i] * x[i] - 1;;
			fx += sin(item);
			*sum += cos(item);
			g[i] += 2 * cos(item)*x[i];
		}
		return fx;
	}
	double EG2
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = 0.5*sin(pow(x[n - 1], 2)) + valgradPoolOneCongestion(EG2Ranged, x, g, n, 0);
		g[n - 1] = cos(pow(x[n - 1], 2))*x[n - 1];
		double item;
		for (int i = block_size - 1; i < n - 1; i += block_size)
		{
			item = x[0] + x[i] * x[i] - 1;;
			fx += sin(item);
			g[0] += cos(item);
			g[i] += 2 * cos(item)*x[i];
		}
		return fx;
	}
	void InitialEG2
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = 0.;
	}
	/*========= 27 ========== ENGVAL1 Function ============= 5000 ==============*/
	double ENGVAL1Ranged
	(
		double *x,
		double *g,
		int begin,
		int end
	)
	{
		double item;
		double fx(0.0);

		for (int i = begin; i < end; i++)
			g[i] = 0.0;

		for (int i = begin; i < end - 1; i++)
		{
			item = x[i] * x[i] + x[i + 1] * x[i + 1];
			fx += item * item + (3 - 4.0*x[i]);
			g[i] += 4.0*item*x[i] - 4.0;
			g[i + 1] += 4.0*item*x[i + 1];
		}
		return fx;
	}
	double ENGVAL1
	(
		double *x,
		double *g,
		const int n
	)
	{
		double item;
		double fx = valgradPool(ENGVAL1Ranged, x, g, n);

		for (int i = block_size - 1; i < n - 1; i += block_size)
		{
			item = x[i] * x[i] + x[i + 1] * x[i + 1];
			fx += item * item + (3 - 4.0*x[i]);
			g[i] += 4.0*item*x[i] - 4.0;
			g[i + 1] += 4.0*item*x[i + 1];
		}
		return fx;
	}
	void InitialENGVAL1
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = 2.0;

	}
	/*===== 28 ======================= EXTROSNB ================== 1000 ================*/
	double EXTROSNBRanged
	(
		double *x,
		double *g,
		int begin,
		int end
	)
	{
		double item;
		double fx = 0;
		for (int i = begin; i < end; i++)
			g[i] = 0.0;
		for (int i = begin + 1; i < end; i++) {
			item = 10 * (x[i] - x[i - 1] * x[i - 1]);
			fx += item * item;
			g[i - 1] += -40.0 * item *x[i - 1];
			g[i] += 20.0 * item;
		}
		return fx;
	}
	double EXTROSNB
	(
		double *x,
		double *g,
		const int n
	)
	{
		double item(x[0] - 1);
		double fx = (pow(item, 2)) + valgradPool(EXTROSNBRanged, x, g, n);
		g[0] += 2 * item;
		for (int i = block_size; i < n; i += block_size) {
			item = 10 * (x[i] - x[i - 1] * x[i - 1]);
			fx += item * item;
			g[i - 1] += -40.0 * item *x[i - 1];
			g[i] += 20.0 * item;
		}
		return fx;
	}
	void InitialEXTROSNB
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = -1.0;
	}
	/*========= 29 ========== FLETCHR Function ============= 1000 ==============*/
	double FLETCHRRanged
	(
		double *x,
		double *g,
		int begin,
		int end
	)
	{
		double item;
		double fx = 0.0;
		for (int i = begin; i < end; i++)
			g[i] = 0;

		for (int i = begin; i < end - 1; i++)
		{
			item = x[i + 1] - x[i] + 1 - x[i] * x[i];
			fx += item * item;
			g[i] += 20.0*item*(-2.0*x[i] - 1.0);
			g[i + 1] += 20.0*item;
		}
		return fx;
	}
	double FLETCHR
	(
		double *x,
		double *g,
		const int n
	)
	{
		double item;
		double fx = valgradPool(FLETCHRRanged, x, g, n);

		for (int i = block_size - 1; i < n - 1; i += block_size)
		{
			item = x[i + 1] - x[i] + 1 - x[i] * x[i];
			fx += item * item;
			g[i] += 20.0*item*(-2.0*x[i] - 1.0);
			g[i + 1] += 20.0*item;
		}
		return 100.0*fx;
	}
	void InitialFLETCHR
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = 0.0;

	}
	/*========= 30 ==============  Freudenstein and Roth: FREUROTH ====== 5000 =======*/
	double FREUROTHRanged
	(
		double *x,
		double *g,
		int begin,
		int end
	)
	{
		double fx = 0.0;
		double item1, item2;
		for (int i = begin; i < end; i++)
			g[i] = 0;
		for (int i = begin; i < end - 1; i++)
		{
			item1 = (-13 + x[i] + ((5 - x[i + 1])*x[i + 1] - 2.0)*x[i + 1]);
			item2 = (-29 + x[i] + ((1 + x[i + 1])*x[i + 1] - 14.0)*x[i + 1]);
			fx += item1 * item1 + item2 * item2;
			g[i] += 2.0*item1 + 2.0*item2;
			g[i + 1] += 2.0*item1*(10 * x[i + 1] - 3.0*x[i + 1] * x[i + 1] - 2.0) +
				2.0*item2*(2 * x[i + 1] + 3.0*x[i + 1] * x[i + 1] - 14.0);
		}
		return fx;
	}
	double FREUROTH
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = valgradPool(FREUROTHRanged, x, g, n);
		double item1, item2;

		for (int i = block_size - 1; i < n - 1; i += block_size)
		{
			item1 = (-13 + x[i] + ((5 - x[i + 1])*x[i + 1] - 2.0)*x[i + 1]);
			item2 = (-29 + x[i] + ((1 + x[i + 1])*x[i + 1] - 14.0)*x[i + 1]);
			fx += item1 * item1 + item2 * item2;
			g[i] += 2.0*item1 + 2.0*item2;
			g[i + 1] += 2.0*item1*(10 * x[i + 1] - 3.0*x[i + 1] * x[i + 1] - 2.0) +
				2.0*item2*(2 * x[i + 1] + 3.0*x[i + 1] * x[i + 1] - 14.0);
		}
		return fx;
	}
	void InitialFREUROTH
	(
		double *x,
		const int n
	)
	{
		x[0] = 0.5;
		x[1] = -2.0;
		for (int i = 2; i < n; i++)
			x[i] = 0.0;
	}

	/*========= 31 ===================  GENHUMPS ====== 5000 ======================*/
	double GENHUMPSRanged
	(
		double *x,
		double *g,
		int begin,
		int end
	)
	{
		for (int i = begin; i < end; i++)
			g[i] = 0;
		double item1 = sin(2.0*x[begin]);
		double item2 = sin(2.0*x[begin + 1]);
		double item11 = item1 * item1;
		double item22 = item2 * item2;
		double t0 = x[begin] * x[begin], t1 = x[begin + 1] * x[begin + 1];
		double fx = item11 * item22 + 0.05*(t0 + t1);
		g[begin] = 4.0*item1*cos(2.0*x[begin]) *item22 + 0.1*x[begin];
		g[begin + 1] = 4.0*item2*cos(2.0*x[begin + 1]) *item11 + 0.1*x[begin + 1];

		for (int i = begin + 1; i < end - 1; i++)
		{
			item1 = item2;
			item2 = sin(2.0*x[i + 1]);
			item11 = item22;
			item22 = item2 * item2;
			t0 = t1;
			t1 = x[i + 1] * x[i + 1];
			fx += item11 * item22 + 0.05*(t0 + t1);
			t0 = 4.0*item1*item2;
			g[i] += t0 * cos(2.0*x[i]) *item2 + 0.1*x[i];
			g[i + 1] += t0 * cos(2.0*x[i + 1])* item1 + 0.1*x[i + 1];
		}
		return fx;
	}
	double GENHUMPS
	(
		double *x,
		double *g,
		const int n
	)
	{
		double item1;
		double item2;
		double item11;
		double item22;
		double t0, t1;
		double fx = valgradPool(GENHUMPSRanged, x, g, n);
		for (int i = block_size - 1; i < n - 1; i += block_size)
		{
			item1 = sin(2.0*x[i]);
			item2 = sin(2.0*x[i + 1]);
			item11 = item1 * item1;
			item22 = item2 * item2;
			t0 = x[i] * x[i];
			t1 = x[i + 1] * x[i + 1];
			fx += item11 * item22 + 0.05*(t0 + t1);
			t0 = 4.0*item1*item2;
			g[i] += t0 * cos(2.0*x[i]) *item2 + 0.1*x[i];
			g[i + 1] += t0 * cos(2.0*x[i + 1])* item1 + 0.1*x[i + 1];
		}
		return fx;
	}
	void InitialGENHUMPS
	(
		double *x,
		const int n
	)
	{
		x[0] = -506.0;
		for (int i = 1; i < n; i++)
			x[i] = 506.2;
	}
	/*===== 32 ======================= GENROSE ====================== 1000 ===================*/
	double GENROSERanged
	(
		double *x,
		double *g,
		int begin,
		int end
	)
	{
		double fx = 0.0;
		double item1, item2;
		for (int i = begin; i < end; i++)
			g[i] = 0;
		for (int i = begin + 1; i < end; i++)
		{
			item1 = 10.0 * (x[i] - x[i - 1] * x[i - 1]);
			item2 = x[i] - 1.0;
			fx += item1 * item1 + item2 * item2;
			g[i - 1] -= 40.0*item1 * x[i - 1];
			g[i] += 20 * item1 + 2 * item2;
		}
		return fx;
	}
	double GENROSE
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = 1.0 + valgradPool(GENROSERanged, x, g, n);
		double item1, item2;

		for (int i = block_size; i < n; i += block_size)
		{
			item1 = 10.0 * (x[i] - x[i - 1] * x[i - 1]);
			item2 = x[i] - 1.0;
			fx += item1 * item1 + item2 * item2;
			g[i - 1] -= 40.0*item1 * x[i - 1];
			g[i] += 20 * item1 + 2 * item2;
		}
		return fx;
	}
	void InitialGENROSE
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = 1.0 / n + 1;
	}
	/*========= 33 ========== LIARWDH Function ============= 5000 ==============*/
	double LIARWDHRanged(double *x, double *g, double *sum, INT index, INT begin, INT end)
	{
		double item1, item2;
		double fx(0.0);
		for (int i = begin; i < end; i++)
		{
			item1 = 2.0*(x[i] * x[i] - x[index]);
			item2 = (x[i] - 1);
			fx += item1 * item1 + item2 * item2;
			g[i] = 8.0*item1*x[i] + 2.0*item2;
			*sum -= 4.0*item1;

		}
		return fx;
	}
	double LIARWDH
	(
		double *x,
		double *g,
		const int n
	)
	{
		return valgradPoolOneCongestion(LIARWDHRanged, x, g, n, 0);
	}
	void InitialLIARWDH
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = 4.0;
	}
	/*========= 34 ========== MOREBV ============= 5000 ==============*/
	double MOREBV
	(
		double *x,
		double *g,
		const int n
	)
	{
		double h(1.0 / n);
		double element, item;
		double fx(0.0);
		for (int i = 0; i < n; i++)
			g[i] = 0;
		//first term
		element = x[1] + h + 1;
		item = 2 * x[1] - x[2] + (h*h*pow(element, 3)) / 2;
		fx += item * item;
		g[1] = 2 * item*(2.0 + h * h*1.5*pow(element, 2));
		g[2] = -2 * item;
		//last term
		element = x[n - 1] + (n - 1)*h + 1;
		item = 2 * x[n - 1] - x[n - 2] + (h*h*pow(element, 3)) / 2;
		fx += item * item;
		g[n - 1] = 4 * item*(2.0 + h * h*1.5*pow(element, 2));
		g[n - 2] = -2 * item;
		for (int i = 2; i < n - 1; i++)
		{
			element = x[i] + i * h + 1;
			item = 2 * x[i] - x[i - 1] - x[i + 1] + (h*h*pow(element, 3)) / 2;
			fx += item * item;
			g[i - 1] -= 2 * item;
			g[i] += 2 * item*(2.0 + h * h*1.5*pow(element, 2));
			g[i + 1] -= 2 * item;
		}
		return fx;
	}
	void InitialMOREBV
	(
		double *x,
		const int n
	)
	{
		double h(1.0 / n);
		x[0] = 0.;
		for (int i = 1; i < n; i++)
			x[i] = i * h*(i*h - 1.0);
	}
	/*===== 35 ======================= NONCVXU2 ============ 1000 ===============*/

	double NONCVXU2
	(
		double *x,
		double *g,
		const int n
	)
	{
		double fx = 0.0;
		for (int i = 0; i < n; i++)
			g[i] = 0;
		double item1, item2;
		for (int j = 0; j < n; j++)
		{
			int i(j + 1);
			item1 = (x[j] + x[(3 * i - 2) % n] + x[(7 * i - 3) % n]);
			item2 = x[j] + x[(3 * i - 2) % n] + x[(7 * i - 3) % n];
			fx += item1 * item1 + 4 * cos(item2);
			g[j] += 2 * item1 - 4 * sin(item2);
			g[(3 * i - 2) % n] += 2 * item1 - 4 * sin(item2);
			g[(7 * i - 3) % n] += 2 * item1 - 4 * sin(item2);
		}

		return fx;
	}
	void InitialNONCVXU2
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = i + 1;
	}
	/*========= 36 ========== NONDIA Function ============= 10000 ==============*/
	double NONDIA
	(
		double *x,
		double *g,
		const int n
	)
	{
		double item(x[0] - 1.0);
		double fx(item*item);

		g[0] = 2.0*item;
		item = 10.0*(x[0] - x[0] * x[0]);
		fx += item * item;
		g[0] += (20.0 - 40.0*x[0])*item;

		for (int i = 1; i < n; i++)
		{
			item = 10.0*(x[0] - x[i] * x[i]);
			fx += item * item;
			g[0] += 20.0*item;
			g[i] = -40.0*x[i] * item;
		}
		return fx;
	}
	void InitialNONDIA
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = -1.0;

	}

	/*========= 37 ========== NONDQUAR Function ============= 10000 ==============*/
	double NONDQUARRanged(double *x, double *g, double *sum, INT index, INT begin, INT end)
	{
		double tmp;
		double fx(0.0);

		for (int i = begin; i < end; i++)
			g[i] = 0.0;

		for (int i = begin; i < end - 2; i++)
		{
			double t0 = x[i] + x[i + 1] + x[index], t1 = t0 * t0, t2 = t1 * t1, t3 = t0 * t1;

			fx += t2;
			tmp = 4.0*t3;
			g[i] += tmp;
			g[i + 1] += tmp;
			*sum += tmp;
		}
		return fx;
	}
	double NONDQUAR
	(
		double *x,
		double *g,
		const int n
	)
	{
		double item, tmp;
		double fx = valgradPoolOneCongestion(NONDQUARRanged, x, g, n, n - 1);

		item = x[0] - x[1];
		fx += item * item;
		g[0] += 2.0 * item;
		g[1] += -2.0 * item;

		item = x[n - 2] + x[n - 1];
		fx += item * item;
		g[n - 2] += 2.0 * item;
		g[n - 1] += 2.0 * item;

		for (int i = block_size - 2; i < n - 2; i += block_size)
		{
			double t0 = x[i] + x[i + 1] + x[n - 1], t1 = t0 * t0, t2 = t1 * t1, t3 = t0 * t1;

			fx += t2;
			tmp = 4.0*t3;
			g[i] += tmp;
			g[i + 1] += tmp;
			g[n - 1] += tmp;
		}
		for (int i = block_size - 1; i < n - 2; i += block_size)
		{
			double t0 = x[i] + x[i + 1] + x[n - 1], t1 = t0 * t0, t2 = t1 * t1, t3 = t0 * t1;

			fx += t2;
			tmp = 4.0*t3;
			g[i] += tmp;
			g[i + 1] += tmp;
			g[n - 1] += tmp;
		}
		return fx;
	}
	void InitialNONDQUAR
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i += 2)
		{
			x[i] = 1.0;
			x[i + 1] = -1.0;
		}

	}
	/*========= 38 ========== PENALTY1 ============= 1000 ==============*/
	double PENALTY1
	(
		double *x,
		double *g,
		const int n
	)
	{
		double item, tmp;
		double tail(0.0);
		double a(1E-5);

		double fx = 0.0;

		for (int i = 0; i < n; i++)
		{
			item = x[i] - 1;
			fx += item * item;
			g[i] = 2.0 * a * item;
			tail += x[i] * x[i];
		}
		fx *= a;

		fx += pow(tail - 0.25, 2);
		tmp = 4 * (tail - 0.25);
		for (int i = 0; i < n; i++)
			g[i] += tmp * x[i];
		return fx;
	}
	void InitialPENALTY1
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = i + 1;

	}
	/*========= 39 ========== PENALTY2 ============= 100 ==============*/
	double PENALTY2
	(
		double *x,
		double *g,
		const int n
	)
	{
		double item1, item2;
		double tail(0.0);
		double a(1E-5);
		double a1(0.2* a);
		const double ExpMinus1by10 = exp(-1 / 10.0);
		for (int i = 0; i < n; i++)
			g[i] = 0.0;
		item1 = x[0] - 0.2;
		double fx(pow(item1, 2));
		g[0] = 2 * item1;
		tail += (n)*pow(x[0], 2);
		double currentTerm = exp(x[0] / 10);
		double prevTer;
		double currentI = exp(1 / 10.0);
		double prevI;
		for (int i = 1; i < n; i++)
		{
			prevTer = currentTerm;
			currentTerm = exp(x[i] / 10);
			prevI = currentI;
			currentI = exp((i + 1) / 10.0);
			item1 = currentTerm + prevTer - currentI - prevI;
			item2 = currentTerm - ExpMinus1by10;
			tail += (n - i)*pow(x[i], 2);
			fx += a * (item1*item1 + item2 * item2);
			g[i - 1] += a1 * item1*prevTer;
			g[i] += a1 * currentTerm*(item1 + item2);
		}
		fx += (tail - 1)*(tail - 1);
		for (int i = 0; i < n; i++)
			g[i] += 4 * (tail - 1)*(n - i)*x[i];
		return fx;
	}

	void InitialPENALTY2
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = 0.5;

	}
	/*========= 40 ========== POWER Function ======== 1000 ==========*/
	double POWERRanged(double *x, double *g, int begin, int end) {
		double item;
		double fx(0.0);

		for (int i = begin; i < end; i++)
		{
			item = (i + 1)*x[i];
			fx += item * item;
			g[i] = 2.0*item*(i + 1);
		}
		return fx;
	}

	double POWER(double *x, double *g, const int n) {
		return valgradPool(POWERRanged, x, g, n);
	}
	void InitialPOWER
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = 1.0;

	}
	/*===== 41 =============== SROSENBR ============= 10000 ================*/
	double SROSENBRRanged
	(
		double *x,
		double *g,
		int begin,
		int end
	)
	{
		int i;
		double fx = 0.0;

		if (begin % 2 != 0) {
			begin++;
		}
		for (i = begin; i < end; i += 2) {
			double t1 = 1.0 - x[i];
			double t2 = 10.0 * (x[i + 1] - x[i] * x[i]);
			g[i + 1] = 20.0 * t2;
			g[i] = -2.0 * (x[i] * g[i + 1] + t1);
			fx += t1 * t1 + t2 * t2;
		}
		return fx;
	}
	double SROSENBR
	(
		double *x,
		double *g,
		const int n
	)
	{
		return valgradPool(SROSENBRRanged, x, g, n);
	}
	void InitialSROSENBR
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i += 2)
		{
			x[i] = -1.2;
			x[i + 1] = 1.0;
		}
	}

	/*========= 42 ========== TRIDIA Function ============= 10000 ==============*/
	double TRIDIARanged
	(
		double *x,
		double *g,
		int begin,
		int end
	)
	{
		double item;
		double fx = 0;
		for (int i = begin; i < end; i++)
			g[i] = 0;

		for (int i = begin + 1; i < end; i++)
		{
			item = (2.0*x[i] - x[i - 1]);
			fx += (i + 1)* item*item;
			g[i] += 4.0*item*(i + 1);
			g[i - 1] -= 2.0*(i + 1)*item;
		}
		return fx;
	}
	double TRIDIA
	(
		double *x,
		double *g,
		const int n
	)
	{
		double item(x[0] - 1.0);
		double fx = (item*item) + valgradPool(TRIDIARanged, x, g, n);
		g[0] += 2.0*item;

		for (int i = block_size; i < n; i += block_size)
		{
			item = (2.0*x[i] - x[i - 1]);
			fx += (i + 1)* item*item;
			g[i] += 4.0*item*(i + 1);
			g[i - 1] -= 2.0*(i + 1)*item;
		}
		return fx;
	}
	void InitialTRIDIA
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = 1.0;
	}
	/*========= 43 =========  Woods ======== 10000 =====================*/
	double WOODSRanged(double *x, double *g, int begin, int end) {
		double fx(0.0);
		double item1, item2, item3, item4, item5, item6;

		while (begin % 4 != 0) {
			begin++;
		}
		for (int i = begin; i < end - 3; i += 4)
		{
			item1 = (x[i + 1] - x[i] * x[i]);
			item2 = (1 - x[i]);
			item3 = (x[i + 3] - x[i + 2] * x[i + 2]);
			item4 = (1 - x[i + 2]);
			item5 = (x[i + 1] + x[i + 3] - 2.0);
			item6 = (x[i + 1] - x[i + 3]);
			fx += (100 * item1*item1 + item2 * item2 + 90 * item3*item3
				+ item4 * item4 + 10.0*(item5*item5 + 0.1*item6*item6));
			g[i] = -400 * item1*x[i] - 2 * item2;
			g[i + 1] = 200 * item1 + 20.0*item5 + 0.2*item6;
			g[i + 2] = -360 * item3*x[i + 2] - 2.0*item4;
			g[i + 3] = 180 * item3 + 20.0*item5 - 0.2*item6;
		}
		return fx;
	}

	double WOODSHelper(double *x, double *g, int n, int i) {
		if (i % 4 != 0) {
			return 0;
		}
		double fx(0.0);
		double item1, item2, item3, item4, item5, item6;
		item1 = (x[i + 1] - x[i] * x[i]);
		item2 = (1 - x[i]);
		item3 = (x[i + 3] - x[i + 2] * x[i + 2]);
		item4 = (1 - x[i + 2]);
		item5 = (x[i + 1] + x[i + 3] - 2.0);
		item6 = (x[i + 1] - x[i + 3]);
		fx += (100 * item1*item1 + item2 * item2 + 90 * item3*item3
			+ item4 * item4 + 10.0*(item5*item5 + 0.1*item6*item6));
		g[i] = -400 * item1*x[i] - 2 * item2;
		g[i + 1] = 200 * item1 + 20.0*item5 + 0.2*item6;
		g[i + 2] = -360 * item3*x[i + 2] - 2.0*item4;
		g[i + 3] = 180 * item3 + 20.0*item5 - 0.2*item6;
		return fx;
	}

	double WOODS(double *x, double *g, const int n) {
		double fx = valgradPool(WOODSRanged, x, g, n);
		for (int i = block_size - 3; i < n - 3; i = i + block_size) {
			fx += WOODSHelper(x, g, n, i);
			fx += WOODSHelper(x, g, n, i + 1);
			fx += WOODSHelper(x, g, n, i + 2);
		}
		return fx;
	}

	void InitialWOODS
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i += 2)
		{
			x[i] = -3;
			x[i + 1] = -1;
		}
	}
	/*========= 40' ========== POWER Function ======== 500  ==========*/
	double POWER1
	(
		double *x,
		double *g,
		const int n
	)
	{
		double item;
		double fx(0.0);

		for (int i = 0; i < n; i += 5)
		{
			item = (i + 1)*x[i];
			fx += item * item;
			g[i] = 2.0*item*(i + 1);

			item = (i + 2)*x[i + 1];
			fx += item * item;
			g[i + 1] = 2.0*item*(i + 2);

			item = (i + 3)*x[i + 2];
			fx += item * item;
			g[i + 2] = 2.0*item*(i + 3);

			item = (i + 4)*x[i + 3];
			fx += item * item;
			g[i + 3] = 2.0*item*(i + 4);

			item = (i + 5)*x[i + 4];
			fx += item * item;
			g[i + 4] = 2.0*item*(i + 5);
		}
		return fx;
	}
	void InitialPOWER1
	(
		double *x,
		const int n
	)
	{
		for (int i = 0; i < n; i++)
			x[i] = 1.0;

	}

	void(*initialize)(double *x, const int n) = InitialPOWER;
	double(*getValGrad)(double *x, double *g, const int n) = POWER;

}