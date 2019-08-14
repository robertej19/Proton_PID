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

public void processFile(String filename,Hist_brbc) {
	HipoDataSource reader = new HipoDataSource()
	reader.open(filename)

	while( reader.hasEvent()){
		DataEvent event = reader.getNextEvent();
		processEvent(event,Hist_brbc)
	}
}

public void processEvent(DataEvent event,Hist_brbc) {
	if(!event.hasBank("REC::Particle")) return
	float startTime = event.getBank("REC::Event").getFloat("startTime",0);
	DataBank reconstructedParticle = event.getBank("REC::Particle");

	p_ind=-1
	if (!hasProton(reconstructedParticle)) return
	if (p_ind>-1){
		(p_momentum, beta_recon,p_theta,p_phi,p_vz,beta_calc,beta_upper,beta_lower) = makeParticle(reconstructedParticle,p_ind)

		DataBank recon_Scint = event.getBank("REC::Scintillator");
		float p_time = 0
		float p_path = 0
		int p_sect = 0
		int p_layer = 0
		//println("Index is: "+recon_Scint.getInt("index",p_ind))
		//println("Detector is: "+recon_Scint.getInt("detector",p_ind))
		if(recon_Scint.getInt("detector",p_ind)==12){
			p_layer = recon_Scint.getInt("layer",p_ind)
			if (![1, 2, 3].contains(p_layer)){
				println("DANGER INVALID LAYER:")
				println(p_layer) }
			p_sect = recon_Scint.getInt("sector",p_ind)
			if (![1, 2, 3, 4, 5, 6].contains(p_sect)){
				println("DANGER INVALID SECTOR:")
				println(p_sect) }
			p_time = recon_Scint.getFloat("time",p_ind)
			p_path = recon_Scint.getFloat("path",p_ind)
			//Question 88 here:
			//fillHists(p_momentum,beta_recon,p_theta,p_phi,p_vz,beta_calc,p_time,p_path)
		}
		else{
			//println("Dectector is not 12, instead it is: "+recon_Scint.getInt("detector",p_ind))
		}

		if (p_momentum < 50){
			if ([1, 2, 3, 4, 5, 6].contains(p_sect)){
				if ([1, 2, 3].contains(p_layer)){
					fillHists(p_momentum,beta_recon,p_theta,p_phi,p_vz,beta_calc,p_time,p_path,p_sect,p_layer,Hist_brbc)
				}
			}
		}

	}
	else return;
}


public boolean hasProton(DataBank reconstructedParticle){
	boolean found = false
	for(int p=0;p<reconstructedParticle.rows();p++){ //Loop over all particles in the event
		if (isProton(reconstructedParticle,p)){ //If we find two or more protons, throw out the event
			//if (found) System.out.println ("Error, two or more Protons found!")
			found=true
		}
	}
	return found
}

public boolean isProton(DataBank reconstructedParticle, int p){
	//if (pID_default_ID_cut(reconstructedParticle,p)&& pID_charge_cut(reconstructedParticle,p)){
	if (pID_beta_momentum_cut(reconstructedParticle,p)&& pID_charge_cut(reconstructedParticle,p)){
	//if (pID_charge_cut(reconstructedParticle,p)){
		p_ind=p //This gives us the index of the row that has the particle event
		return true
	}
	else return false
}

public boolean pID_default_ID_cut(DataBank reconstructedParticle, int p){
  if(reconstructedParticle.getInt("pid",p)==2212) return true;
  else return false;
}

public boolean pID_beta_momentum_cut(DataBank reconstructedParticle, int p){
	(p_momentum, beta_recon,p_theta,p_phi,p_vz,beta_calc,beta_upper,beta_lower) = makeParticle(reconstructedParticle,p)

	if((beta_recon<beta_upper) && (beta_recon>beta_lower)){
		 //println("Momentum value is: "+p_momentum)
		 //println("Beta Recon value is: "+beta_recon)
		 //println("Beta Calc value is: "+beta_calc)
		 return true
	 }
  else return false;

}

public boolean pID_charge_cut(DataBank reconstructedParticle, int p){
  if(reconstructedParticle.getInt("charge",p)==1) return true;
  else return false;
}

def makeParticle(DataBank reconstructedParticle,int p_ind){
		//println("p_ind p_ind is: "+p_ind)
		//"REC::Scintillator"
		float px = reconstructedParticle.getFloat("px",p_ind)
		float py = reconstructedParticle.getFloat("py",p_ind)
		float pz = reconstructedParticle.getFloat("pz",p_ind)
		float beta_recon = reconstructedParticle.getFloat("beta",p_ind)
		float p_momentum = (float)Math.sqrt(px*px+py*py+pz*pz)
		float p_vz = reconstructedParticle.getFloat("vz",p_ind)
		float p_vx = reconstructedParticle.getFloat("vx",p_ind)
		float p_vy = reconstructedParticle.getFloat("vy",p_ind)
		Ve = new LorentzVector(px,py,pz,p_momentum)
		float p_phi = (float) Math.toDegrees(Ve.phi())
		float p_theta = (float) Math.toDegrees(Ve.theta())
		float p_mass = 0.938 //Proton mass in GeV

		float scale_factor = 0.30
		float beta_calc = (float)Math.sqrt(p_momentum*p_momentum/(p_momentum*p_momentum+p_mass*p_mass))
		float p_mom_up = p_momentum*(1+scale_factor)
		float p_mom_low = p_momentum*(1-scale_factor)
		float beta_upper = (float)Math.sqrt(p_mom_up*p_mom_up/(p_mom_up*p_mom_up+p_mass*p_mass))
		float beta_lower = (float)Math.sqrt(p_mom_low*p_mom_low/(p_mom_low*p_mom_low+p_mass*p_mass))

		return [p_momentum, beta_recon,p_theta,p_phi,p_vz,beta_calc,beta_upper,beta_lower]
}

public void fillHists(p_momentum,beta_recon,p_theta,p_phi,p_vz,beta_calc,p_time,p_path,p_sect,p_layer,Histogram_List){
	//int hist_layer = 6*(p_layer-1)+p_sect-1
	"""
	int hist_layer = p_sect-1
	H_proton_beta_momentum[hist_layer].fill(p_momentum,beta_recon)
	H_proton_DeltaBeta_momentum[hist_layer].fill(p_momentum,beta_recon-beta_calc)
	H_proton_mom[hist_layer].fill(p_momentum);
	H_beta_recon_beta_calc[hist_layer].fill(beta_recon-beta_calc);
	H_proton_vz_mom[hist_layer].fill(p_momentum,p_vz);
	H_proton_theta_mom[hist_layer].fill(p_momentum,p_theta)
	H_proton_phi_mom[hist_layer].fill(p_momentum,p_phi)
	H_proton_theta_phi[hist_layer].fill(p_phi,p_theta);

	H_proton_time[hist_layer].fill(p_time);
	H_proton_path[hist_layer].fill(p_path);
	//H_proton_sect.fill(p_sect);

	//for(int isec=1;isec<=6;isec++)
	//for(int ilay=1;ilay<=3;ilay++)
	//println "Trying to fill histogram"""
	Histogram_List["sec${p_sect}_layer${p_layer}"].fill(beta_recon-beta_calc);
	//Hist_brbc["sec2_layer2"].fill(1);


}

"""------------------------ Variable Definitions -------------------------"""

int run = args[0].split("/")[-1].split('\\.')[0][-4..-1].toInteger()

float EB = 10.6f
if(run>6607) EB=10.2f

int p_ind, p_sect

TDirectory out = new TDirectory()
out.mkdir('/'+run)
out.cd('/'+run)

def Hist_momentum 			= [:].withDefault{new H1F("hist_momentum${it}"				, "Momentum ${it}"										,100,0,EB)}
def Hist_time 					= [:].withDefault{new H1F("hist_time${it}"						, "Time ${it}"												,100,0,250)}
def Hist_path_length 		= [:].withDefault{new H1F("Hist_path_length${it}"			, "Path Length ${it}"									,100,400,1000)}
def Hist_vz 						= [:].withDefault{new H1F("hist_vz${it}"							, "Z vertex ${it}"										,100,-25,25)}
def Hist_beta_recon			= [:].withDefault{new H1F("hist_beta_recon${it}"			, "REC::Part Beta vs Calc Beta ${it}"	,100,-1,1)}
def Hist_beta_p 				= [:].withDefault{new H2F("hist_beta_p${it}"					, "Beta vs. Momentum ${it}"						,100,0,EB,100,0,1)}
def Hist_deltaB_p 			= [:].withDefault{new H2F("Hist_deltaB_p${it}"				, "Delta B vs. Momentum ${it}"				,400,0,EB,400,-0.2,0.2)}
def Hist_momentum_vz 		= [:].withDefault{new H2F("hist_momentum_vz${it}"			, "Momentum vs. Vz ${it}"							,100,0,EB,100,-25,25)}
def Hist_momentum_theta = [:].withDefault{new H2F("hist_momentum_theta${it}"	, "Momentum vs. Theta ${it}"					,100,0,EB,100,0,40)}
def Hist_momentum_phi 	= [:].withDefault{new H2F("hist_momentum_phi${it}"		, "Momentum vs. Phi ${it}"						,100,0,EB,100,-180, 180)}
def Hist_theta_phi 			= [:].withDefault{new H2F("hist_theta_phi${it}"				, "Theta vs. Phi ${it}"								,100,-180, 180,100,0,40)}


Hist_List = [Hist_momentum,Hist_time,Hist_path_length,Hist_vz,Hist_beta_recon,Hist_beta_p]

println Hist_List[1]

"""------------------------ Start of Program -------------------------"""

filenum=-1 //There should be able to get rid of this filenum issue
for (arg in args){
	filenum=filenum+1
	if (filenum==0) continue
	processFile(arg,Hist_beta_recon)
}

"""
(0..<max_hists).each{
	out.addDataSet(H_proton_beta_momentum[it])
	out.addDataSet(H_proton_DeltaBeta_momentum[it])
	out.addDataSet(H_proton_mom[it])
	out.addDataSet(H_beta_recon_beta_calc[it])
	out.addDataSet(H_proton_vz_mom[it])
	out.addDataSet(H_proton_theta_mom[it])
	out.addDataSet(H_proton_phi_mom[it])
	out.addDataSet(H_proton_theta_phi[it])
	out.addDataSet(H_proton_time[it])
	out.addDataSet(H_proton_path[it])
}"""

for(int isec=1;isec<=6;isec++){
 for(int ilay=1;ilay<=3;ilay++){
   out.addDataSet(Hist_beta_recon["sec${isec}_layer${ilay}"])
 }
}

out.writeFile('proton_pID_new'+run+'.hipo')
