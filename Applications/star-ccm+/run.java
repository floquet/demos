// STAR-CCM+ macro: run.java
package macro;

import java.util.*;

import star.common.*;
import star.base.neo.*;


public class run extends StarMacro {

  public void execute() {

    Simulation simulation_0 = 
      getActiveSimulation();

    StepStoppingCriterion stepStoppingCriterion_0 = 
      ((StepStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Steps"));

    stepStoppingCriterion_0.setMaximumNumberSteps(2);

    simulation_0.getSimulationIterator().run();

  }
}
