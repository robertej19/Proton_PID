import java.io.*;
import java.util.*;
import org.jlab.groot.data.TDirectory
import org.jlab.groot.data.GraphErrors
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.math.F1D;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.clas.physics.Vector3;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.graphics.EmbeddedCanvas;

"""------------------------ Function Definitions -------------------------"""

public void fillHists(p_momentum,p_theta){
	H_proton_theta_momentum[p_sect-1].fill(p_momentum,p_theta)
}

public void processFile(String filename) {
	HipoDataSource reader = new HipoDataSource()
	reader.open(filename)

	while( reader.hasEvent()){
		DataEvent event = reader.getNextEvent();
		processEvent(event)
	}
}
public void processEvent(DataEvent event) {
	if(!event.hasBank("REC::Particle")) return
	DataBank particleBank = event.getBank("REC::Particle");
	e_index=-1
	if (!hasElectron(particleBank)) return
	if (e_index>-1){
		p_momentum, p_theta = makeElectron(particleBank)
		fillHists(p_momentum,p_theta)
	}
	else return;
}


public boolean hasElectron(DataBank recPart){
	boolean found = false
	for(int p=0;p<recPart.rows();p++){
		if (isElectron(recPart,p)){
			if (found) System.out.println ("Error, two or more electrons found!")
			found=true
		}
	}
	return found
}
public boolean isElectron(DataBank recPart, int p){
	if (pID_default_electron_cut(recPart,p)&& pID_charge_cut(recPart,p)){
		e_index=p
		return true
	}
	else return false
}

public boolean pID_default_electron_cut(DataBank recPart, int p){
  if(recPart.getInt("pid",p)==2212) return true;
  else return false;
}
public boolean pID_charge_cut(DataBank recPart, int p){
  if(recPart.getInt("charge",p)==1) return true;
  else return false;
}

public void makeElectron(DataBank recPart){
		int ei=e_index
		println("e_index ei is: "+ei)
		float px = recPart.getFloat("px",ei)
		float py = recPart.getFloat("py",ei)
		float pz = recPart.getFloat("pz",ei)
		float p_momentum = (float)Math.sqrt(px*px+py*py+pz*pz)
		e_vz = recPart.getFloat("vz",ei)
		e_vx = recPart.getFloat("vx",ei)
		e_vy = recPart.getFloat("vy",ei)
		e_phi = (float)Math.toDegrees(Math.atan2(py,px))
		if(e_phi<0) e_phi+=360;
		e_phi=360-e_phi;
		e_phi=e_phi-150;
		if (e_phi<0) e_phi+=360;
		p_sect = (int) Math.ceil(e_phi/60);
		Ve = new LorentzVector(px,py,pz,p_momentum)
		e_phi = (float) Math.toDegrees(Ve.phi())
		float p_theta = (float) Math.toDegrees(Ve.theta())
		return p_momentum, p_theta
}

"""------------------------ Variable Definitions -------------------------"""

def run = args[0].toInteger()
float EB = 10.6f
if(run>6607) EB=10.2f

int e_index, p_sect
float e_phi, e_vx, e_vy, e_vz
LorentzVector Ve = new LorentzVector()

TDirectory out = new TDirectory()
out.mkdir('/'+run)
out.cd('/'+run)

H_proton_theta_momentum =(0..<6).collect{
	def h1 = new H2F("H_proton_theta_momentum_S"+(it+1), "H_proton_theta_momentum_S"+(it+1),100,0,EB,100,0,40);
	return h1
}

"""------------------------ Start of Program -------------------------"""

filenum=-1 //There should be able to get rid of this filenum issue
for (arg in args){
	filenum=filenum+1
	if (filenum==0) continue
	processFile(arg)
}

(0..<6).each{
	out.addDataSet(H_proton_theta_momentum[it])
}

out.writeFile('proton_pID_new'+run+'.hipo')
