void Main::enzo_finalize (Simulation * simulation)
{
  const Config * config   = simulation->config();
  Monitor * monitor = simulation->monitor();

  simulation->finalize();
  
  int cycle_final = config->testing_cycle_final;

  unit_class ("Enzo-E");
  unit_func  ("final cycle");
  if (cycle_final != 0) {
    unit_assert (simulation->cycle()==cycle_final);
    monitor->print ("Testing","actual   cycle:  %d",simulation->cycle());
    monitor->print ("Testing","expected cycle:  %d",cycle_final);
  }

  double time_tolerance = config->testing_time_tolerance;

  unit_class ("Enzo-E");
  unit_func  ("final time");
  monitor->print ("Testing","actual   time:  %.15g",simulation->time());
  monitor->print ("Testing","tolerance:      %g",time_tolerance);
  
  if (config->testing_time_final.size() > 0 &&
      config->testing_time_final[0] > 0.0) {
    double err_rel_min = std::numeric_limits<double>::max();
    for (size_t i=0; i<config->testing_time_final.size(); i++) {
      double time_final=config->testing_time_final[i];
      double err_rel = cello::err_rel(simulation->time(),time_final);
      err_rel_min = std::min(err_rel_min,err_rel);
      monitor->print ("Testing","expected time:  %.15g",time_final);
      monitor->print ("Testing","relative error: %g",err_rel);
    }
    monitor->print ("Testing","minimum relative error: %g",err_rel_min);
    unit_assert ( err_rel_min < time_tolerance);
  }

  monitor->print ("","END ENZO-E");

}

