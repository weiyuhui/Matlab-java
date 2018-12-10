package pyulear;

public class test_float {
	
	private static final float PI = 3.1415926f;

	public static void main(String[] args) throws Exception {
		
		float[] data = new float[500];
		float fs = 1;
		float f0 = 0.3F;
		for(int i=0;i<data.length;i++){
			data[i] = (float) Math.sin(2*PI*f0*1/fs*i);
			//data[i] = (float)i;
			//System.out.println(data[i]);
		}
		
		//len 指数据的长度
		int len = 256;

		
		Pyulear_float pyulear_ = new Pyulear_float(data, len);
		float[] redata_ = Pyulear_float.freqMe(data, len, true);
		System.out.println(redata_.length);
		for(int j =0;j<redata_.length;j++){
			System.out.println(redata_[j]);
		}
		
	}

}
