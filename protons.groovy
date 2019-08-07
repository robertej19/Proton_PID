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
	p_ind=-1
	if (!hasProton(particleBank)) return
	if (p_ind>-1){
		(p_momentum, beta_groot) = makeProton(particleBank,p_ind)
		float tofTime =	event.getBank("FTOF::hits").getFloat("time",p_ind)
		//float pl =	event.getBank("FTOF::hits").getFloat("path",p_ind)
		//println("Pathlength is: "+pl)
		//float tof = tofTime - startTime
		//float dis = 2
		//float beta = dis/tof
		//println("Beta: "+beta)
		fillHists(p_momentum,beta_groot)
	}
	else return;
}


public boolean hasProton(DataBank reconstructedParticle){
	boolean found = false
	for(int p=0;p<reconstructedParticle.rows();p++){
		if (isProton(reconstructedParticle,p)){
			//if (found) System.out.println ("Error, two or more Protons found!")
			found=true
		}
	}
	return found
}
public boolean isProton(DataBank reconstructedParticle, int p){
	if (pID_default_ID_cut(reconstructedParticle,p)&& pID_charge_cut(reconstructedParticle,p)){
		p_ind=p
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

public void fillHists(p_momentum,beta){
	H_proton_beta_momentum[p_sect-1].fill(p_momentum,beta)
	H_proton_mom[p_sect-1].fill(p_momentum);
	H_proton_vz_mom[p_sect-1].fill(p_momentum,p_vz);
	H_proton_theta_mom[p_sect-1].fill(p_momentum,p_theta)
	H_proton_phi_mom[p_sect-1].fill(p_momentum,p_phi)
	H_proton_theta_phi[p_sect-1].fill(p_phi,p_theta);
}


float p_phi, p_vx, p_vy, p_vz
LorentzVector Ve = new LorentzVector()

def makeProton(DataBank reconstructedParticle,int p_ind){
		//println("p_ind p_ind is: "+p_ind)
		float px = reconstructedParticle.getFloat("px",p_ind)
		float py = reconstructedParticle.getFloat("py",p_ind)
		float pz = reconstructedParticle.getFloat("pz",p_ind)
		float beta = reconstructedParticle.getFloat("beta",p_ind)
		println("Groot beta: "+beta)
		float p_momentum = (float)Math.sqrt(px*px+py*py+pz*pz)
		p_vz = reconstructedParticle.getFloat("vz",p_ind)
		p_vx = reconstructedParticle.getFloat("vx",p_ind)
		p_vy = reconstructedParticle.getFloat("vy",p_ind)
		p_phi = (float)Math.toDegrees(Math.atan2(py,px))
		if(p_phi<0) p_phi+=360;
		p_phi=360-p_phi;
		p_phi=p_phi-150;
		if (p_phi<0) p_phi+=360;
		p_sect = (int) Math.ceil(p_phi/60);
		Ve = new LorentzVector(px,py,pz,p_momentum)
		p_phi = (float) Math.toDegrees(Ve.phi())
		float p_theta = (float) Math.toDegrees(Ve.theta())
		float p_mass = 0.938 //Proton mass in GeV
		float rbeta = (float)Math.sqrt(p_momentum*p_momentum/p_mass/p_mass/(1+p_momentum*p_momentum/p_mass/p_mass))
		return [p_momentum, rbeta]
}

"""------------------------ Variable Definitions -------------------------"""

def run = args[0].toInteger()
float EB = 10.6f
if(run>6607) EB=10.2f

int p_ind, p_sect

TDirectory out = new TDirectory()
out.mkdir('/'+run)
out.cd('/'+run)

H_proton_beta_momentum =(0..<6).collect{
	def h1 = new H2F("H_proton_beta_momentum_S"+(it+1), "H_proton_beta_momentum_S"+(it+1),300,0,EB,100,0,1);
	return h1
}

H_proton_mom =(0..5).collect{
	def h1 = new H1F("H_proton_mom_S"+(it+1), "H_proton_mom_S"+(it+1),100, 0, EB);
	return h1
}

H_proton_vz_mom =(0..<6).collect{
	def h1 = new H2F("H_proton_vz_mom_S"+(it+1), "H_proton_vz_mom_S"+(it+1),100,0,EB,100,-25,25);
	return h1
}

H_proton_theta_mom =(0..<6).collect{
	def h1 = new H2F("H_proton_theta_mom_S"+(it+1), "H_proton_theta_mom_S"+(it+1),100,0,EB,100,0,40);
	return h1
}

H_proton_phi_mom = (0..<6).collect{
	def h1 = new H2F("H_proton_phi_mom_S"+(it+1), "H_proton_phi_mom_S"+(it+1),100,0,EB,100,-180,180);
	return h1
}

H_proton_theta_phi =(0..<6).collect{
	def h1 = new H2F("H_proton_theta_phi_S"+(it+1), "H_proton_theta_phi_S"+(it+1),100,-180,180,100,0,40);
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
	out.addDataSet(H_proton_mom[it])
	out.addDataSet(H_proton_vz_mom[it])
	out.addDataSet(H_proton_theta_mom[it])
	out.addDataSet(H_proton_phi_mom[it])
	out.addDataSet(H_proton_theta_phi[it])
}

out.writeFile('proton_pID_new'+run+'.hipo')
