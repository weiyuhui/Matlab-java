package test;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;


import interpolation.Interpolation;
import interpolation.PointStruct;
public class Testspline {
// 根据X排序
class SortByX implements Comparator {
	public int compare(Object obj1, Object obj2) {
		PointStruct point1 = (PointStruct) obj1;
		PointStruct point2 = (PointStruct) obj2;
		if (point1.dx > point2.dx)
                 return 1;
             else
                 return 0;
         }
     } 
     //转化成数组
     public static ArrayList<Double> getArrayfList(List<PointStruct> l, String type) {
         ArrayList<Double> tmplist = new ArrayList<Double>();
         for (PointStruct p : l) {
             if (type == "ix") {
                 tmplist.add((double) p.ix);
             } else if (type == "iy") {
                 tmplist.add((double) p.iy);
             } else if (type == "dx") {
                 tmplist.add(p.dx);
             } else if (type == "dy") {
                 tmplist.add(p.dy);
             }
         }
         return tmplist;
     }
     
     public static void main(String[] args) {
         
         //将来做插值的点列表
         List <PointStruct> pl=new  ArrayList<PointStruct>();
         
         //【1】 模拟 sin函数值  做初始化特定点（插值点或样条），  曲线将来就穿过下面几点    
         //        0.1    1    2    3    4         //        0.099833417    0.841470985    0.909297427    0.141120008    0.61802
         
         pl.add(new PointStruct(0.1,0.099833417));
         pl.add(new PointStruct(1,0.841470985));
         pl.add(new PointStruct(2,0.909297427));
         pl.add(new PointStruct(3,0.141120008));
         pl.add(new PointStruct(4,0.61802));
      
         //【2】 添加 需要的未知点（新插值点）， 为了使曲线光滑  采用小的步长 填充在已知点之间
         List <PointStruct> target=new  ArrayList<PointStruct>();        
         double step=0.1;        
         for(int k=0;k<pl.size();k++){
             
             if((k+1)<pl.size()){
                 double tmpd=0;
                 while(tmpd<pl.get(k+1).dx){
                     tmpd=tmpd+step;
                     target.add(new PointStruct(tmpd,0.0));
                 }
             }
         }
         
         //把点集合转化成为 x,y 各自的坐标点 Double[]集合 
         ArrayList tmp_x = getArrayfList(pl,"dx");
         ArrayList tmp_y = getArrayfList(pl,"dy");
         Double[] xspline=new Double[tmp_x.size()];
         Double[] yspline=new Double[tmp_x.size()];        
         for(int e=0;e<tmp_x.size();e++){
             xspline[e]=(Double) tmp_x.get(e);
             yspline[e]=(Double) tmp_y.get(e);
         } 
         Double[] x = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
         Double[] y = {0.0,0.841470984807897,0.909297426825682,0.141120008059867,-0.756802495307928,-0.958924274663139,-0.279415498198926,0.656986598718789,0.989358246623382,0.412118485241757,-0.544021110889370};
         double[] xx = {0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,};
         
         //【3】根据 插值点初始化插值对象,并用spline样条函数进行计算（做三次样条插值运算）
         Interpolation ip = new Interpolation(x, y);
         for(int j=0;j<xx.length;j++){
             double yy =ip.spline(xx[j]);
             System.out.println(yy);
         }
         
//         //把新生成的样条添加到点集合中  做以后画图用 
//         pl.addAll(target);         
//         //根据x点的大小进行点集合的排序
//         Testspline t2=new Testspline();
//         Collections.sort(pl, t2.new SortByX());
//         
//         //因为结果太小，所有必要放大一定度数，同时 Java画图最像是整数，所有把双精度转成整形
//         int intScale=100; 
//         for(PointStruct p:pl){
//              p.ix=(int) (p.dx*intScale);
//              p.iy=(int) (p.dy*intScale);
//         }
//         
//         //查看中间结果 即画图所用到的所有点
//         int times=0;
//         for(PointStruct p:pl){
//             System.out.println(" order： :"+times++ +"  p.ix:"+p.ix+"    p.iy:"+p.iy);
//         }
//         
////         //【4】 调用画图类 画图
//         //http://www.cnblogs.com/rojas/p/4595509.html
////         Drawlineforspline dlfs=new Drawlineforspline(pl); 
////         Mypanel myp=dlfs.new Mypanel();
////         dlfs.add(myp); 
 
     }
 }