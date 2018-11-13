package test;

public class testJama {
	
//	public static void main(String[] args) {
//		//double[][] arr1 = {{-1, 1, 0}, {-4, 3, 0}, {1, 0, 2}};
//		double[][] arr1 = {{-1, 1, 0}};
//		double[][] arr2 = {{1, 1, 0}, {1, 3, 0}, {1, 0, 2}};
//		DoubleMatrix A = new DoubleMatrix(arr1);
//		DoubleMatrix B = new DoubleMatrix(arr2);
//		DoubleMatrix F = B.divColumnVector(A);
//		F.print();
//		
//	}
	
	public static void main(String[] args) {
	double[] x = { 1,8,15,23,31,39,47,55,63,71,79,87,94,102,110,118,126,134,142,149,157,165,173,181,188,196,204,212,220,227,235 };
	double[] y = {1.34026,1.3402,1.2658,1.2500,1.2822,1.2658,1.2658,1.2500,1.2822,1.2658,1.2658,1.2658,1.2658,1.28205,1.2658,1.2658,1.2987,1.2658,1.2658,1.2658,1.2822,1.2822,1.2987,1.2822,1.2987,1.2500,1.2821,1.2821,1.2821,1.2821,1.2658};
	double[] xx = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
	
//	double[] x = { 1,8,15,23,31,39};
//	double[] y = {1.34026,1.3402,1.2658,1.2500,1.2822,1.2658};

	try {
		double[] v = Spline_Mat.spline(x, y, xx);
		
		for(int i = 0; i<v.length; i++){
			System.out.println(v[i]);
		}
	} catch (Exception e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}


	}

}
