import org.jlab.groot.data.TDirectory
import org.jlab.groot.data.GraphErrors
import org.jlab.groot.data.H1F
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;


for(arg in args) {
  TDirectory dir = new TDirectory()
  dir.readFile(arg)
  println("Printing directory")
  dir.pwd()
  dir.getDirectoryList()
  println(dir.pwd())
  dir.ls()
  def h1 = dir.getObject(String.format('/5038/H_proton_beta_momentum_S3'))
  //...
//
//  def h_n = dir.getObject(String.format('/path/to/histogram_n'))

  EmbeddedCanvas canvas  = new EmbeddedCanvas();
  canvas.setSize(2500,1500);
  //canvas.divide(sqrt(n),sqrt(n));
  canvas.setAxisTitleSize(40);
  canvas.setAxisFontSize(40);
  canvas.setTitleSize(40);
  //canvas.cd(0);
	canvas.draw(h1);
  //...
  //canvas.cd(n-1);canvas.draw(h_n);
  canvas.save('dumped_histogram'+arg+'.png')
}
