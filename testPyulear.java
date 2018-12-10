package pyulear;

public class testPyulear {

	private static final double PI = 3.1415926;

	public static void main(String[] args) throws Exception {
		
//		double[] x = { 1, 2, 3, 4, 5, 6, 7 };
//		double[] data = {-0.00129471153178508, -0.782583239790595, 0.719703362593308, -0.239663085896841,
//				0.0239142302628919, -0.119973473168448, 0.208764401325136, -0.0243503692351991, -0.121620813721807,
//				-0.354759387832616, 1.13916150561757, -1.21851284327035, 0.404726126155552, 0.465056513324389,
//				-0.745718136433046, 0.440972026096388, 0.151840157370897, -0.564604371731284, 0.451174017125206,
//				-0.0422719500206666, -0.0997973074083268, -0.254813358866707, 0.766184971338706, -0.895762676355103,
//				0.508612842210098, 0.0676286113547881, -0.558694114421898, 0.861563457559211, -0.760078898190947,
//				0.156453222381764, 0.467179811236297, -0.483225394948457, 0.0135903649044826, 0.283809780953458,
//				-0.295484603229648, 0.572921858148884, -1.21093490144773, 1.32161369110346, -0.272618408007197,
//				-1.06504462600085, 1.38952846749448, -0.721542451432558, -0.0197767812857384, 0.145331904171351,
//				0.547520518528617, -1.38917263428019, 1.35576249846226, -0.489459881696096, -0.435831234953442,
//				0.876363975503758, -0.460212859300983, -0.390745661887404, 0.458153955804443, 0.175358957759666,
//				-0.190788121714871, -0.330373212340608, 0.380082020537307, -0.198774022077930, 0.231270884731384,
//				0.0312308562055841, -0.653696027893262, 0.881174385373959, -0.475389257568191, -0.0475815696674845,
//				0.309562015994878, -0.387104248413835, 0.419740464632498, -0.364141732727299, 0.150279296063090,
//				0.0395264189716402, 0.110185380850692, -0.450915238300992, 0.504853962143243, -0.312270802078884,
//				0.114429634330873, 0.0779212060475069, 0.292407279326990, -1.43106188381539, 1.69107520675785,
//				-0.00371435307263707 };
		
		double[] data = new double[300];
		double fs = 10;
		double f0 = 0.4F;
		for(int i=0;i<data.length;i++){
			data[i] = Math.sin(2*PI*f0*1/fs*i);
			//data[i] = (double)i;
			//System.out.println(data[i]);
		}
		
		//len ָ���ݵĳ���
		int len = 512;
//		Pyulear_8 pyulear = new Pyulear_8(data, len);
//		double[] redata = Pyulear_8.freqMe(data, len, true);
//		System.out.println(redata.length);
//		for(int j =0;j<redata.length;j++){
//			System.out.println(redata[j]);
//		}
		
		//��C++����У��֮��ARģ��ϵ����ͬ��
		PyulearBurgC pyulear_ = new PyulearBurgC(data, len);
		double[] redata_ = PyulearBurgC.freqMe(data, len, true);
		//System.out.println(redata_.length);
		for(int j =0;j<redata_.length;j++){
			System.out.println(redata_[j]);
		}
	
		//��Դ����վ�ϵ�C���ԣ���ʦ������ַ
//		Pyulear_C pyulear_ = new Pyulear_C(data, len);
//		double[] redata_ = pyulear_.freqMe(data, len, true);
//		//System.out.println(redata_.length);
//		for(int j =0;j<redata_.length;j++){
//			System.out.println(redata_[j]);
//		}
		
	}

}