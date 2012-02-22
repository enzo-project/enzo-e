#ifdef CONFIG_USE_CHARM
void Main::enzo_finalize (Simulation * simulation)
#else
void enzo_finalize (Simulation * simulation)
#endif
{

  Parameters * parameters = simulation->parameters();
  Monitor * monitor       = simulation->monitor();

  simulation->finalize();
  
  // parameter: Testing : cycle_final

  int    cycle_final = parameters->value_integer("Testing:cycle_final",0);

  unit_class ("Enzo-P");
  unit_func  ("final cycle");
  if (cycle_final != 0) {
    unit_assert (simulation->cycle()==cycle_final);
    monitor->print ("Testing","actual   cycle:  %d",simulation->cycle());
    monitor->print ("Testing","expected cycle:  %d",cycle_final);
  }

  // parameter: Testing : time_final

  double time_final  = parameters->value_float("Testing:time_final",0.0);

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
}

