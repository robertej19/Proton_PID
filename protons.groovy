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

public void fillHists(p_momentum,beta){
	H_proton_beta_momentum[p_sect-1].fill(p_momentum,beta)
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
	float startTime = event.getBank("REC::Event").getFloat("startTime",0);
	DataBank particleBank = event.getBank("REC::Particle");
	//DataBank forwardTOFBank = event.getBank("REC::Particle")
	//if(event.hasBank("REC::Event"))STT =
	//println(starttime)
	e_index=-1
	if (!hasElectron(particleBank)) return
	if (e_index>-1){
		(p_momentum, beta_groot) = makeElectron(particleBank,e_index)
		float tofTime =	event.getBank("FTOF::hits").getFloat("time",e_index)
		//float pl =	event.getBank("FTOF::hits").getFloat("path",e_index)
		//println("Pathlength is: "+pl)
		float tof = tofTime - startTime
		float dis = 2
		float beta = dis/tof
		//println("Beta: "+beta)
		fillHists(p_momentum,beta_groot)
	}
	else return;
}


public boolean hasElectron(DataBank reconstructedParticle){
	boolean found = false
	for(int p=0;p<reconstructedParticle.rows();p++){
		if (isElectron(reconstructedParticle,p)){
			//if (found) System.out.println ("Error, two or more electrons found!")
			found=true
		}
	}
	return found
}
public boolean isElectron(DataBank reconstructedParticle, int p){
	if (pID_default_ID_cut(reconstructedParticle,p)&& pID_charge_cut(reconstructedParticle,p)){
		e_index=p
		return true
	}
	else return false
}

public boolean pID_default_ID_cut(DataBank reconstructedParticle, int p){
  if(reconstructedParticle.getInt("pid",p)==2212) return true;
  else return false;
}
public boolean pID_charge_cut(DataBank reconstructedParticle, int p){
  if(reconstructedParticle.getInt("charge",p)==1) return true;
  else return false;
}

float e_phi, e_vx, e_vy, e_vz
LorentzVector Ve = new LorentzVector()

def makeElectron(DataBank reconstructedParticle,int e_index){
		int ei=e_index
		//println("e_index ei is: "+ei)
		float px = reconstructedParticle.getFloat("px",ei)
		float py = reconstructedParticle.getFloat("py",ei)
		float pz = reconstructedParticle.getFloat("pz",ei)
		float beta = reconstructedParticle.getFloat("beta",ei)
		println("Groot beta: "+beta)
		float p_momentum = (float)Math.sqrt(px*px+py*py+pz*pz)
		e_vz = reconstructedParticle.getFloat("vz",ei)
		e_vx = reconstructedParticle.getFloat("vx",ei)
		e_vy = reconstructedParticle.getFloat("vy",ei)
		e_phi = (float)Math.toDegrees(Math.atan2(py,px))
		if(e_phi<0) e_phi+=360;
		e_phi=360-e_phi;
		e_phi=e_phi-150;
		if (e_phi<0) e_phi+=360;
		p_sect = (int) Math.ceil(e_phi/60);
		Ve = new LorentzVector(px,py,pz,p_momentum)
		e_phi = (float) Math.toDegrees(Ve.phi())
		float p_theta = (float) Math.toDegrees(Ve.theta())
		float p_mass = 0.938 //Proton mass in GeV
		float rbeta = (float)Math.sqrt(p_momentum*p_momentum/p_mass/p_mass/(1+p_momentum*p_momentum/p_mass/p_mass))
		return [p_momentum, rbeta]
}

"""------------------------ Variable Definitions -------------------------"""

def run = args[0].toInteger()
float EB = 10.6f
if(run>6607) EB=10.2f

int e_index, p_sect

TDirectory out = new TDirectory()
out.mkdir('/'+run)
out.cd('/'+run)

H_proton_beta_momentum =(0..<6).collect{
	def h1 = new H2F("H_proton_beta_momentum_S"+(it+1), "H_proton_beta_momentum_S"+(it+1),300,0,EB,100,0,1);
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
	out.addDataSet(H_proton_beta_momentum[it])
}

out.writeFile('proton_pID_new'+run+'.hipo')
