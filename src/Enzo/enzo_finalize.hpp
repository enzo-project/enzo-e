void Main::enzo_finalize (Simulation * simulation)
{
  const Config * config   = simulation->config();
  Monitor * monitor = simulation->monitor();

  simulation->finalize();
  
  int cycle_final = config->testing_cycle_final;

  unit_class ("Enzo-P");
  unit_func  ("final cycle");
  if (cycle_final != 0) {
    unit_assert (simulation->cycle()==cycle_final);
    monitor->print ("Testing","actual   cycle:  %d",simulation->cycle());
    monitor->print ("Testing","expected cycle:  %d",cycle_final);
  }

  double time_final  = config->testing_time_final;

  unit_class ("Enzo-P");
  unit_func  ("final time");
  if (time_final != 0.0) {
    double err_rel = cello::err_rel(simulation->time(),time_final);
    double mach_eps = cello::machine_epsilon(precision_default);
    unit_assert ( err_rel < 100*mach_eps);
    monitor->print ("Testing","actual   time:  %.15g",simulation->time());
    monitor->print ("Testing","expected time:  %.15g",time_final);
    monitor->print ("Testing","relative error: %g",err_rel);
    monitor->print ("Testing","100*mach_eps:   %g",100*mach_eps);
  }

  monitor->print ("","END ENZO-P");

}

