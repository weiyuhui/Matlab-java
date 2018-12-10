package pyulear;

public class PyulearBurgC {

	private double[] data;
	private static int rst_len;
	private static int AR_LEAST_ORDER = 60;

	public PyulearBurgC(double[] data, int rst_len) {
		this.data = data;
		this.rst_len = rst_len;
	}

	public static double[] freqMe(double[] data, int rst_len, boolean bRemoveDC) {
		// TODO Auto-generated method stub
		double[] prstdata = new double[rst_len];
		if (data == null | data.length < AR_LEAST_ORDER) {
			return null;
		}
		double[] fpe = new double[AR_LEAST_ORDER];
		// 初始化值---
		double pm = 0;
		double[] g = new double[AR_LEAST_ORDER];
		int k = 0;
		// AR_LEAST_ORDER原始值是多少？----------60,循环遍历，找到最合适的阶数
		int lg = AR_LEAST_ORDER - 10;
		double[] pdata_buf = null;
		if (bRemoveDC == true) {
			pdata_buf = new double[data.length];
			double dc = 0;
			for (int j = 0; j < data.length; j++) {
				dc += data[j];
			}
			dc = dc / data.length;
			for (int j = 0; j < data.length; j++) {
				pdata_buf[j] = data[j] - dc;
			}
		}


		if (bRemoveDC == true) {
			data = pdata_buf;
		}
		/*
		 * burg 过程 因为函数有多个输出，所以运用过程来代替函数
		 */
		int ic = 1;
		int im = 1;
		int i, j, jj, n, nn, jn1, mpnn = 0, mmnn = 0, jnn1, nnj1, lg1, jk1, lext;
		double dm, sum, ftemp, sn, sd, t1, t2;
		double[] am, pef, per;
		double[] h = new double[AR_LEAST_ORDER];

		int skip_flag = 0;
		int m = data.length - 1;
		// if (m < 5)
		// return 0;

		am = new double[m * 2 + 200];
		pef = new double[m * 2 + 200];
		per = new double[m * 2 + 200];

		lext = lg + 1;
		sum = 0.;
		for (i = 1; i <= m; i++) {
			sum = sum + data[i] * data[i];
		}
		// m = data_len -1;返回数据的长度
		am[1] = sum / m;
		pm = am[1];
		dm = pm;
		if (im == 0) {
			if (ic == 2)
				fpe[1] = (double) (m + 2) / (double) (m - 2) * pm;
			else
				fpe[1] = (double) (m + 1) / (double) (m - 1) * pm;
		} else {
			if (ic == 2)
				fpe[1] = (double) (m + 4) / (double) (m - 4) * pm;
			else
				fpe[1] = (double) (m + 2) / (double) (m - 2) * pm;
		}
		ftemp = fpe[1];
		fpe[1] = 0.;

		for (nn = 2; nn <= lg; nn++) {
			skip_flag = 0;

			n = nn - 2;
			if (n == 0) {
				for (j = 1; j <= m; j++) {
					pef[j] = 0.0;
					per[j] = 0.0;
				}
			}
			sn = 0.;
			sd = 0.;
			jj = m - n - 1;

			for (j = 1; j <= jj; j++) {
				jn1 = j + n + 1;
				t1 = data[jn1];
				t2 = data[j];
				sn = sn - (2.0 * (t1 + pef[j]) * (t2 + per[j]));
				sd = sd + (Math.pow((t1 + pef[j]), 2.) + Math.pow((t2 + per[j]), 2.));
			}
			if (sd == 0.f) {
				am = null;
				pef = null;
				per = null;
				// return 0;
			}
			g[nn] = sn / sd;
			if (Math.abs(g[nn]) > 1.) // goto aj2;
			{
				skip_flag = 2;
				break;
			}
			if (n == 0) // goto aj3;
			{
				for (j = 1; j <= jj; j++) {
					jnn1 = j + nn - 1;
					per[j] = per[j] + g[nn] * pef[j] + g[nn] * data[jnn1];
					pef[j] = pef[j + 1] + g[nn] * per[j + 1] + g[nn] * data[j + 1];
				}
				sum = 0.;
				for (j = 2; j <= nn; j++) {
					nnj1 = nn + 1 - j;
					sum = sum - am[nnj1] * g[j];
				}
				am[nn] = sum;
				dm = (float) ((1. - Math.pow(g[nn], 2)) * dm);
				pm = dm;
				if (nn == m)
					continue;

				skip_flag = 0;
				if (ic != 2) // goto aj5;
				{
					mpnn = m + nn;
					mmnn = m - nn;
					if (im != 0) // goto aj6;
					{
						mpnn++;
						mmnn--;
					}
					// goto aj6;
					skip_flag = 6;
				}

				if (skip_flag != 6) {
					mpnn = m + nn * 2;
					mmnn = m - nn * 2;
					if (im != 0) // goto aj6;
					{
						mpnn = mpnn + 2;
						mmnn = mmnn - 2;
					}
				}

				if (mmnn == 0) {
					am = null;
					pef = null;
					per = null;
					// return 0;
				}
				// 精度计算有误差：0.458289 0.452812
				fpe[nn] = (double) mpnn / (double) mmnn * pm;
				if (ftemp == 0.f) {
					am = null;
					pef = null;
					per = null;
					// return 0;
				}
				fpe[nn] = fpe[nn] / ftemp;
				fpe[nn] = Math.log10(fpe[nn]);

				skip_flag = 0;
			} else {

				// n==0;
				for (j = 1; j <= n; j++) {
					k = n - j + 2;
					h[j + 1] = g[j + 1] + g[nn] * g[k];
				}
				for (j = 1; j <= n; j++) {
					g[j + 1] = h[j + 1];
				}
				jj--;

				for (j = 1; j <= jj; j++) {
					jnn1 = j + nn - 1;
					per[j] = per[j] + g[nn] * pef[j] + g[nn] * data[jnn1];
					pef[j] = pef[j + 1] + g[nn] * per[j + 1] + g[nn] * data[j + 1];
				}
				sum = 0.;
				for (j = 2; j <= nn; j++) {
					nnj1 = nn + 1 - j;
					sum = sum - am[nnj1] * g[j];
				}
				am[nn] = sum;
				dm = (float) ((1. - Math.pow(g[nn], 2)) * dm);
				pm = dm;
				if (nn == m)
					continue;

				skip_flag = 0;
				if (ic != 2) // goto aj5;
				{
					mpnn = m + nn;
					mmnn = m - nn;
					if (im != 0) // goto aj6;
					{
						mpnn++;
						mmnn--;
					}
					// goto aj6;
					skip_flag = 6;
				}

				if (skip_flag != 6) {
					mpnn = m + nn * 2;
					mmnn = m - nn * 2;
					if (im != 0) // goto aj6;
					{
						mpnn = mpnn + 2;
						mmnn = mmnn - 2;
					}
				}

				if (mmnn == 0) {
					am = null;
					pef = null;
					per = null;
					// return 0;
				}
				// 精度计算有误差：0.458289 0.452812
				fpe[nn] = (double) mpnn / (double) mmnn * pm;
				if (ftemp == 0.f) {
					am = null;
					pef = null;
					per = null;
					// return 0;
				}
				fpe[nn] = fpe[nn] / ftemp;
				fpe[nn] = Math.log10(fpe[nn]);

				skip_flag = 0;
			}
		}

		if (skip_flag != 2) {
			g[1] = 1.0;
			lg1 = lg + 1;
			for (j = lg1; j <= lext; j++) {
				sum = 0.;
				for (k = 2; k <= lg; k++) {
					jk1 = j + 1 - k;
					sum = sum - am[jk1] * g[k];
				}
				am[j] = sum;
			}
			am = null;
			pef = null;
			per = null;
			// return 1;
		}

		// 确定了最后的阶数
		// 2.burg

		double fpemin = fpe[2];
		int naropt = 2;
		for (int j2 = 2; j2 <= lg; j2++) {// for 50 times, AR order
			if (fpe[j2] <= fpemin) {
				fpemin = fpe[j2];
				naropt = j2 - 1;
			}
		}

		if (bRemoveDC == true) {
			data = pdata_buf;
		}

		lg = naropt + 1;
		am = null;
		pef = null;
		per = null;
		pm = 0;
		fpe = new double[AR_LEAST_ORDER];
		// 初始化值---
		pm = 0;
		g = new double[AR_LEAST_ORDER];
		// return 2;
		ic = 1;
		im = 1;
		mpnn = 0;
		mmnn = 0;
		// double dm, sum, ftemp, sn, sd, t1, t2;
		// am, pef, per;
		h = new double[AR_LEAST_ORDER];

		skip_flag = 0;
		m = data.length - 1;
		// if (m < 5)
		// return 0;

		am = new double[m * 2 + 200];
		pef = new double[m * 2 + 200];
		per = new double[m * 2 + 200];

		lext = lg + 1;
		sum = 0.;
		for (i = 1; i <= m; i++) {
			sum = sum + data[i] * data[i];
		}
		// m = data_len -1;返回数据的长度
		am[1] = sum / m;
		pm = am[1];
		dm = pm;
		if (im == 0) {
			if (ic == 2)
				fpe[1] = (double) (m + 2) / (double) (m - 2) * pm;
			else
				fpe[1] = (double) (m + 1) / (double) (m - 1) * pm;
		} else {
			if (ic == 2)
				fpe[1] = (double) (m + 4) / (double) (m - 4) * pm;
			else
				fpe[1] = (double) (m + 2) / (double) (m - 2) * pm;
		}
		ftemp = fpe[1];
		fpe[1] = 0.;

		for (nn = 2; nn <= lg; nn++) {
			skip_flag = 0;

			n = nn - 2;
			if (n == 0) {
				for (j = 1; j <= m; j++) {
					pef[j] = 0.0;
					per[j] = 0.0;
				}
			}
			sn = 0.;
			sd = 0.;
			jj = m - n - 1;

			for (j = 1; j <= jj; j++) {
				jn1 = j + n + 1;
				t1 = data[jn1];
				t2 = data[j];
				sn = sn - (2.0 * (t1 + pef[j]) * (t2 + per[j]));
				sd = sd + (Math.pow((t1 + pef[j]), 2.) + Math.pow((t2 + per[j]), 2.));
			}
			if (sd == 0.f) {
				am = null;
				pef = null;
				per = null;
				// return 0;
			}
			g[nn] = sn / sd;
			if (Math.abs(g[nn]) > 1.) // goto aj2;
			{
				skip_flag = 2;
				break;
			}
			if (n == 0) // goto aj3;
			{
				for (j = 1; j <= jj; j++) {
					jnn1 = j + nn - 1;
					per[j] = per[j] + g[nn] * pef[j] + g[nn] * data[jnn1];
					pef[j] = pef[j + 1] + g[nn] * per[j + 1] + g[nn] * data[j + 1];
				}
				sum = 0.;
				for (j = 2; j <= nn; j++) {
					nnj1 = nn + 1 - j;
					sum = sum - am[nnj1] * g[j];
				}
				am[nn] = sum;
				dm = (float) ((1. - Math.pow(g[nn], 2)) * dm);
				pm = dm;
				if (nn == m)
					continue;

				skip_flag = 0;
				if (ic != 2) // goto aj5;
				{
					mpnn = m + nn;
					mmnn = m - nn;
					if (im != 0) // goto aj6;
					{
						mpnn++;
						mmnn--;
					}
					// goto aj6;
					skip_flag = 6;
				}

				if (skip_flag != 6) {
					mpnn = m + nn * 2;
					mmnn = m - nn * 2;
					if (im != 0) // goto aj6;
					{
						mpnn = mpnn + 2;
						mmnn = mmnn - 2;
					}
				}

				if (mmnn == 0) {
					am = null;
					pef = null;
					per = null;
					// return 0;
				}
				// 精度计算有误差：0.458289 0.452812
				fpe[nn] = (double) mpnn / (double) mmnn * pm;
				if (ftemp == 0.f) {
					am = null;
					pef = null;
					per = null;
					// return 0;
				}
				fpe[nn] = fpe[nn] / ftemp;
				fpe[nn] = Math.log10(fpe[nn]);

				skip_flag = 0;
			} else {

				// n==0;
				for (j = 1; j <= n; j++) {
					k = n - j + 2;
					h[j + 1] = g[j + 1] + g[nn] * g[k];
				}
				for (j = 1; j <= n; j++) {
					g[j + 1] = h[j + 1];
				}
				jj--;

				for (j = 1; j <= jj; j++) {
					jnn1 = j + nn - 1;
					per[j] = per[j] + g[nn] * pef[j] + g[nn] * data[jnn1];
					pef[j] = pef[j + 1] + g[nn] * per[j + 1] + g[nn] * data[j + 1];
				}
				sum = 0.;
				for (j = 2; j <= nn; j++) {
					nnj1 = nn + 1 - j;
					sum = sum - am[nnj1] * g[j];
				}
				am[nn] = sum;
				dm = (float) ((1. - Math.pow(g[nn], 2)) * dm);
				pm = dm;
				if (nn == m)
					continue;

				skip_flag = 0;
				if (ic != 2) // goto aj5;
				{
					mpnn = m + nn;
					mmnn = m - nn;
					if (im != 0) // goto aj6;
					{
						mpnn++;
						mmnn--;
					}
					// goto aj6;
					skip_flag = 6;
				}

				if (skip_flag != 6) {
					mpnn = m + nn * 2;
					mmnn = m - nn * 2;
					if (im != 0) // goto aj6;
					{
						mpnn = mpnn + 2;
						mmnn = mmnn - 2;
					}
				}

				if (mmnn == 0) {
					am = null;
					pef = null;
					per = null;
					// return 0;
				}
				// 精度计算有误差：0.458289 0.452812
				fpe[nn] = (double) mpnn / (double) mmnn * pm;
				if (ftemp == 0.f) {
					am = null;
					pef = null;
					per = null;
					// return 0;
				}
				fpe[nn] = fpe[nn] / ftemp;
				fpe[nn] = Math.log10(fpe[nn]);

				skip_flag = 0;
			}
		}

		if (skip_flag != 2) {
			g[1] = 1.0;
			lg1 = lg + 1;
			for (j = lg1; j <= lext; j++) {
				sum = 0.;
				for (k = 2; k <= lg; k++) {
					jk1 = j + 1 - k;
					sum = sum - am[jk1] * g[k];
				}
				am[j] = sum;
			}
			am = null;
			pef = null;
			per = null;
			// return 1;
		}
		// return 2;

		double[] me_para = new double[AR_LEAST_ORDER];

		me_para[0] = pm;
		for (k = 1; k <= naropt; k++) {
			me_para[k] = g[k + 1];
		}

		/*
		 * 
		 * 1.看看Matlab的频段划分是如何实现的 2.缩减选择阶数的那个环节； 3.整理返回数据的长度与AR阶数的关系
		 * A:0-0.5之间选择频率成分，等分成len（预定义）段，然后根据数据长度的个数，选择频段。。
		 */
		prstdata = CalData(me_para, naropt, 0., 0.5);

		double sum1, sum2;

		if (bRemoveDC == true) {
			sum1 = 0.f;
			for (int j3 = 0; j3 < data.length; j3++) {
				sum1 += pdata_buf[j3] * pdata_buf[j3];
			}
		} else {
			sum1 = 0;
			for (int j3 = 0; j3 < data.length; j3++) {
				sum1 += data[j3] * data[j3];
			}
		}

		sum2 = 0;
		for (int j3 = 0; j3 < rst_len; j3++) {
			sum2 += prstdata[j3];
		}

		double coeff;
		if (sum2 != 0.f)
			coeff = sum1 / sum2;
		else
			coeff = 1.f;

		for (int j3 = 0; j3 < rst_len; j3++) {
			prstdata[j3] = prstdata[j3] * coeff;
		}

		// if (pdata_buf!=NULL)
		// delete []pdata_buf;

		return prstdata;
	}

	private static double[] CalData(double[] porgpara, int paralen, double f0, double f1) {
		// TODO Auto-generated method stub
		if (f0 < 0.)
			f0 = 0.;
		if (f0 > 0.49)
			f0 = 0.49;
		if (f1 < 0.01)
			f1 = 0.01;
		if (f1 > 0.5)
			f1 = 0.5;
		if (f1 < f0)
			f1 = f0 + 0.01;

		int datalen = rst_len;
		// int paralen = porgpara.length;//返回数值的长度，自己定义
		double[] pdata = new double[datalen];

		double f = f0;
		double det_f = ((f1 - f0) / datalen);
		int nf = datalen;
		double sumr, sumi, tmpf;
		int j, k;
		for (j = 1; j < nf; j++) {
			tmpf = (2 * Math.PI * f);
			sumr = 1.;
			sumi = 0.;
			for (k = 1; k <= paralen; k++) {
				sumr = sumr + (porgpara[k] * Math.cos(tmpf * k));
				sumi = sumi - (porgpara[k] * Math.sin(tmpf * k));
			}
			tmpf = (sumr * sumr + sumi * sumi);
			if (tmpf > 0)
				pdata[j - 1] = (porgpara[0] / (sumr * sumr + sumi * sumi));
			else
				pdata[j - 1] = 0;
			f += det_f;
		}

		return pdata;

	}

}
