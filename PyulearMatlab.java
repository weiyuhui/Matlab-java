package pyulear;

public class PyulearMatlab {
	
	public static void main(String[] args) {
		
		double[] data = new double[300];
		double fs = 10;
		double f0 = 0.3F;
		for(int i=0;i<data.length;i++){
			data[i] = Math.sin(2*Math.PI*f0*1/fs*i);
			//data[i] = (double)i;
			//System.out.println(data[i]);
		}
		
		//len 指数据的长度
		int len = 512;
		double[] filterData = freqMe(data, len, true);
		for(int j =0;j<filterData.length;j++){
			System.out.println(filterData[j]);
		}
		
		
	}
	

	private static double[] freqMe(double[] data, int len, boolean b) {
		// TODO Auto-generated method stub
		int order = 8;
		boolean isReal = true;
		
		return null;
	}

}
