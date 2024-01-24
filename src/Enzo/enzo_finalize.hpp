void Main::enzo_finalize (Simulation * simulation)
{
  const Config * config   = simulation->config();
  Monitor * monitor = simulation->monitor();

  simulation->finalize();
  
  int cycle_final = config->testing_cycle_final;

  unit_class ("Enzo-E");
  unit_func  ("final cycle");
  if (cycle_final != 0) {
    const int cycle = simulation->state()->cycle();
    unit_assert (cycle==cycle_final);
    monitor->print ("Testing","actual   cycle:  %d",cycle);
    monitor->print ("Testing","expected cycle:  %d",cycle_final);
  }

  double time_tolerance = config->testing_time_tolerance;

  unit_class ("Enzo-E");
  unit_func  ("final time");
  const double time = simulation->state()->time();
  monitor->print ("Testing","actual sim-time:  %.15g", time);
  monitor->print ("Testing","tolerance:      %g",time_tolerance);
  
  if (config->testing_time_final.size() > 0 &&
      config->testing_time_final[0] > 0.0) {
    double err_rel_min = std::numeric_limits<double>::max();
    for (size_t i=0; i<config->testing_time_final.size(); i++) {
      double time_final=config->testing_time_final[i];
      double err_rel = cello::err_rel(time,time_final);
      err_rel_min = std::min(err_rel_min,err_rel);
      monitor->print ("Testing","expected sim-time:  %.15g",time_final);
      monitor->print ("Testing","relative error: %g",err_rel);
    }
    monitor->print ("Testing","minimum relative error: %g",err_rel_min);
    unit_assert ( err_rel_min < time_tolerance);
  }

  monitor->print ("","END ENZO-E");

}

