#!/usr/bin/groovy

import java.io.*
import java.util.*
import org.jlab.groot.data.TDirectory
import org.jlab.groot.data.GraphErrors
import org.jlab.groot.group.DataGroup
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.math.F1D
import org.jlab.groot.fitter.DataFitter
import org.jlab.io.base.DataBank
import org.jlab.io.base.DataEvent
import org.jlab.io.hipo.HipoDataSource
import org.jlab.io.hipo.HipoDataSync
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Vector3
import org.jlab.clas.physics.LorentzVector
import org.jlab.groot.base.GStyle
import org.jlab.groot.graphics.EmbeddedCanvas

float EB = 10.6f
int run = args[0].split("/")[-1].split('\\.')[0][-4..-1].toInteger()
if(run>6607) EB=10.2f

def reader = new HipoDataSource()
reader.open(args[0])

def Hist_momentum 			= [:].withDefault{new H1F("Hist_momentum${it}"				, "Momentum ${it}"										,100,0,EB)}
def Hist_time 					= [:].withDefault{new H1F("Hist_time${it}"						, "Time ${it}"												,100,0,250)}
def Hist_path_length 		= [:].withDefault{new H1F("Hist_path_length${it}"			, "Path Length ${it}"									,100,400,1000)}
def Hist_vz 						= [:].withDefault{new H1F("Hist_vz${it}"							, "Z vertex ${it}"										,100,-25,25)}
def Hist_beta_recon			= [:].withDefault{new H1F("Hist_beta_recon${it}"			, "REC::Part Beta vs Calc Beta ${it}"	,100,-1,1)}
def Hist_beta_p 				= [:].withDefault{new H2F("Hist_beta_p${it}"					, "Beta vs. Momentum ${it}"						,100,0,EB,100,0,1)}
def Hist_deltaB_p 			= [:].withDefault{new H2F("Hist_deltaB_p${it}"				, "Delta B vs. Momentum ${it}"				,400,0,EB,400,-0.2,0.2)}
def Hist_momentum_vz 		= [:].withDefault{new H2F("Hist_momentum_vz${it}"			, "Momentum vs. Vz ${it}"							,100,0,EB,100,-25,25)}
def Hist_momentum_theta = [:].withDefault{new H2F("Hist_momentum_theta${it}"	, "Momentum vs. Theta ${it}"					,100,0,EB,100,0,40)}
def Hist_momentum_phi 	= [:].withDefault{new H2F("Hist_momentum_phi${it}"		, "Momentum vs. Phi ${it}"						,100,0,EB,100,-180, 180)}
def Hist_theta_phi 			= [:].withDefault{new H2F("Hist_theta_phi${it}"				, "Theta vs. Phi ${it}"								,100,-180, 180,100,0,40)}

while(reader.hasEvent()) {
//for (int i=0; i < 15; i++) {
	//println("event i is: $i")
  def event = reader.getNextEvent()
	if(!event.hasBank("REC::Particle")) continue

	//float startTime = event.getBank("REC::Event").getFloat("startTime",0);
	DataBank reconstructedParticle = event.getBank("REC::Particle")

	for(int p_ind=0;p_ind<event.getBank("REC::Particle").rows();p_ind++){ //Loop over all particles in the event
		//println("P is $p_ind")
		if(!reconstructedParticle.getInt("charge",p_ind)==1) continue
		//println("positive charge, continuing")
		//if(!reconstructedParticle.getInt("pid",p)==2212) return false;
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

		//if(beta_recon<beta_upper || beta_recon>beta_lower) return

		DataBank recon_Scint = event.getBank("REC::Scintillator")
		if(recon_Scint.getInt("detector",p_ind)==12){
			p_layer = recon_Scint.getInt("layer",p_ind)
			p_sect = recon_Scint.getInt("sector",p_ind)
			p_time = recon_Scint.getFloat("time",p_ind)
			p_path = recon_Scint.getFloat("path",p_ind)

			if ([1, 2, 3, 4, 5, 6].contains(p_sect) && [1, 2, 3].contains(p_layer)){
				title = "sec${p_sect}_layer${p_layer}"
				Hist_momentum[title].fill(p_momentum)
				Hist_time[title].fill(p_time)
				Hist_path_length[title].fill(p_path)
				Hist_vz[title].fill(p_vz)
				Hist_beta_recon[title].fill(beta_recon-beta_calc)
				Hist_beta_p[title].fill(p_momentum,beta_recon)
				Hist_deltaB_p[title].fill(p_momentum,beta_recon-beta_calc)
				Hist_momentum_vz[title].fill(p_momentum,p_vz)
				Hist_momentum_theta[title].fill(p_momentum,p_theta)
				Hist_momentum_phi[title].fill(p_momentum,p_phi)
				Hist_theta_phi[title].fill(p_phi,p_theta)
			}
		}
		else{
			println("Dectector is not 12, instead it is: "+recon_Scint.getInt("detector",p_ind))
		}
	}
}

reader.close()

TDirectory out = new TDirectory()
out.mkdir('/'+run)
out.cd('/'+run)

for(int isec=1;isec<=6;isec++){
 for(int ilay=1;ilay<=3;ilay++){
	 title = "sec${isec}_layer${ilay}"
	 out.addDataSet(Hist_momentum[title])
	 out.addDataSet(Hist_time[title])
	 out.addDataSet(Hist_path_length[title])
	 out.addDataSet(Hist_vz[title])
	 out.addDataSet(Hist_beta_recon[title])
	 out.addDataSet(Hist_beta_p[title])
	 out.addDataSet(Hist_deltaB_p[title])
	 out.addDataSet(Hist_momentum_vz[title])
	 out.addDataSet(Hist_momentum_theta[title])
	 out.addDataSet(Hist_momentum_phi[title])
	 out.addDataSet(Hist_theta_phi[title])
 }
}

out.writeFile('proton_pID_new'+run+'.hipo')
