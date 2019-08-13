#!/usr/bin/groovy

import org.jlab.io.hipo.HipoDataSource
import org.jlab.io.hipo.HipoDataSync
import org.jlab.detector.base.DetectorType

def reader = new HipoDataSource()
reader.open(args[0])


for (int i=0; i < 5; i++) {
  def event = reader.getNextEvent()
  if (event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter")) {
    def event_cal = event.getBank("REC::Calorimeter")
    def pind_array = event_cal.getShort('pindex')*.toInteger()
    def sect_array = event_cal.getByte('sector')
    println "pindex array for event $i is: $pind_array and has sector array: $sect_array"

  }
}

reader.close()
