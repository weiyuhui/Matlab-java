package filtfilt;

import java.util.ArrayList;

public class FilterFilterArray {
	
	
	/**
	 * 
	 * @param B
	 *            滤波器参数B
	 * @param A
	 *            滤波器参数A
	 * @param input
	 *            要滤波的数组 ： 要求是一维，长度要大于B数组。
	 * @return 滤波后的结果 返回double数组。
	 */
	public static double[] doFiltfilt(double[] B,double[] A, double[] X) {

		int len = X.length;
		int na = A.length;
		int nb = B.length;
		int nfilt = (nb > na) ? nb : na;
		int nfact = 3 * (nfilt - 1);
		if (len <= nfact)
			throw new RuntimeException("输入数值X长度太小，数据最少是滤波器阶数的三倍");
		resize(B,nfilt, 0);
		resize(A,nfilt, 0);
		
		int[] rows = null;
		int[] cols = null;

		add_index_range(rows, 0, nfilt - 2,1);
		if (nfilt > 2)
		{
			add_index_range(rows, 1, nfilt - 2, 1);
			add_index_range(rows, 0, nfilt - 3, 1);
		}
		add_index_const(cols, 0, nfilt - 1);
		if (nfilt > 2)
		{
			add_index_range(cols, 1, nfilt - 2,1);
			add_index_range(cols, 1, nfilt - 2,1);
		}
		int klen = rows.length;

		double[] data = new double[klen];
		data[0] =  1 + A[1];
		int j = 1;
		if (nfilt > 2)
		{
			for (int i = 2; i < nfilt; i++)
				data[j++] = A[i];
			for (int i = 0; i < nfilt - 2; i++)
				data[j++] = 1.0;
			for (int i = 0; i < nfilt - 2; i++)
				data[j++] = -1.0;
		}
		double[] leftpad = subvector_reverse(X, nfact, 1);
		changeArray2(leftpad,2*X[0]);
		
		double[] rightpad = subvector_reverse(X, len - 2, len - nfact - 1);
		changeArray2(rightpad,2*X[len - 1]);
		
		double y0;
		double[] signal1 = null;
		double[] signal2 = null;
		double[] zi = null;
		//reserve(signal1,leftpad.size() + X.size() + rightpad.size(),0);
		append_vector(signal1, leftpad);
		append_vector(signal1, X);
		append_vector(signal1, rightpad);
		
		double [][] sp = Zeros(max_val(rows)+1,max_val(cols)+1); 
		for (int k = 0; k < klen; ++k)
			sp[rows[k]][cols[k]] = data[k];
		double[]bb = B;
		double[]aa = A;
		//Mat.inv(sp)
		//bb.segment(1, nfilt - 1) - (bb(0) * aa.segment(1, nfilt - 1))
		
		double[][] ZZi = Mat.multi(Mat.inv(sp),calc(segment(bb,1,nfilt - 1),bb[0],segment(aa,1,nfilt - 1)));
		
		resize(zi,ZZi.length,1);
		
		changeZi(ZZi,zi,signal1[0]);
		filter(B, A, signal1, signal2, zi);
		reverse(signal2);
		changeZi(ZZi,zi,signal2[0]);
		filter(B, A, signal2, signal1, zi);
		double[] Y = subvector_reverse(signal1, signal1.length - nfact - 1, nfact);
		return Y;

	}

	private static void reverse(double[] signal2) {
		int i=0;
		int j=signal2.length-1;
		while(i<j){
			swap(signal2,i,j);
			i++;
			j--;
		}
	}

	private static void swap(double[] signal2, int i, int j) {
		double temp = signal2[j];
		signal2[j] = signal2[i];
		signal2[i] = temp;
	}

	private static void changeZi(double[][] zZi, double[] zi,
			Double double1) {
		for (int i = 0; i < zZi.length; i++) {
			zi[i] = zZi[i][0]*double1;
		}
		
	}

	private static double[][] calc(double[] segment, double d, double[] segment2) {
		double[][] ret = new double[segment.length][1];
		for (int i = 0; i < segment.length; i++) {
			ret[i][0] = segment[i]-d*segment2[i];
		}
		return ret;
	}

	private static double[] segment(double[] bb, int i, int j) {
		double[]ret=new double[j-i+1];
		for (int k = 0; k < j-i+1; k++) {
			ret[k] = bb[i+k];
		}
		return ret;
	}

//	private static double[] map(ArrayList<Double> b) {
//		double[] ret = new double[b.size()];
//		for (int i = 0; i < ret.length; i++) {
//			ret[i] = b.get(i);
//		}
//		return ret;
//	}

	private static double[][] Zeros(int ii, int jj) {
		double [][] sp = new double[ii][jj];
		for (int i = 0; i < ii; i++)
			for (int j = 0; j < jj; j++)
				sp[i][j] = 0;
		return sp;
	}

	public static void filter(double[] B, double[] A,double[] X, double[] Y,double[] Zi) {
		if (A.length == 0)
			throw new RuntimeException("A 数组为空！");
		boolean flagA = true;
		for (Double doubleA : A) {
			if (doubleA != 0) {
				flagA = false;
			}
		}
		if (flagA) {
			throw new RuntimeException("A 数组至少要有一个数不为零！");
		}
		if (A[0] == 0) {
			throw new RuntimeException("A 数组第一个元素不能为零！");
		}
		changeArray(A, A[0]);
		changeArray(B, A[0]);
		
		
		int input_size = X.length;
		int filter_order = max(A.length,B.length);
		resize(B,filter_order,0);
		resize(A,filter_order,0);
		resize(Zi,filter_order,0);
		resize(Y,input_size,0);
		
		for (int i = 0; i < input_size; i++) {
			int order = filter_order - 1;
			while(order!=0){
				if(i>=order)
					Zi[order-1] = B[order]*X[i-order]-A[order]*Y[i-order]+Zi[order];
				--order;
			}
			Y[i] = B[0]*X[i]+Zi[0];
		}		
		Zi = remove(Zi);		
	}

	



	private static double[] remove(double[] d) {
		double[] re = new double[d.length-1];
		for (int i = 0; i < re.length; i++)
			re[i] = d[i];
		// TODO Auto-generated method stub
		return re;
	}
	
	private static double[] add(double[] d, double data) {
		if(d == null){
			double[] da = {data};
			return da;
		}else{
			double[] re = new double[d.length+1];
			for (int i = 0; i < re.length; i++)
				re[i] = d[i];
			re[d.length] = data;
			return re;	
		}
		
	}
	
	private static int[] add(int[] d, int data) {
		if(d == null){
			int[] da = {data};
			return da;
		}else{
			int[] re = new int[d.length+1];
			for (int i = 0; i < re.length; i++)
				re[i] = d[i];
			re[d.length] = data;
			return re;
		}
		
	}

	private static double[] resize(double[] a, int i, double j) {
		if(a.length>=i){
			return a;
		}else{
			int size = a.length;
			for (int j2 = size; j2 < i; j2++) {
				a = add(a, j);
			}
			return a;
		}
	}
	
	

	private static int max(int size, int size2) {
		if(size >size2)
			return size;
		else
			return size2;
	}

	static void changeArray(double[] vec, double a0) {
		for (int i = 0; i < vec.length; i++) {
			vec[i] = vec[i]/a0;
		}
	}
	
	static void changeArray2(double[] vec, double a0) {
		for (int i = 0; i < vec.length; i++) {
			vec[i] = a0-vec[i];
		}
	}

	static void add_index_range(int[] indices, int beg, int end, int inc) {
		if(indices == null){
			int[] data = new int[end-beg+1];
			for(int i=0;i<data.length;i++){
				data[i] = i;
			}
			indices = data;
		}else{
			for(int i = beg; i <= end; i += inc){
				indices = add(indices, i);
			}
		}
	}

	static void add_index_const(int[] indices, int value, int numel) {
		//if(indices == null){
			int[] da = new int[numel];
			while (numel-- != 0)
				da[numel] = value;
			indices = da;
		//}
		
	}

	static void append_vector(double[] vec, double[] tail) {
		if(vec == null){
			vec = tail;
		}else{
			for(int i=0;i<tail.length;i++){
				vec = add(vec, tail[i]);
			}
		}
	}

	static double[] subvector_reverse(double[] vec, int idx_end, int idx_start) {
		double[] resultArrayList = new double[idx_end- idx_start + 1];
		for (int i = 0; i < idx_end- idx_start + 1; i++) {
			resultArrayList[i] = 0.0;
		}
		int endindex = idx_end - idx_start;
		for (int i = idx_start; i <= idx_end; i++)
			resultArrayList[endindex--] = vec[i];
		return resultArrayList;
	}

	static int max_val(int[] vec) {
		int temp = vec[0];
		for (Integer integer : vec) {
			if (temp < integer)
				temp = integer;
		}
		return temp;
	}


}
