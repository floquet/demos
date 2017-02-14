// STAR-CCM+ macro: star_dem.java
package macro;

import java.util.*;

import star.common.*;
import star.base.neo.*;


public class star_dem extends StarMacro {

  public void execute() {

    Simulation simulation_0 = 
      getActiveSimulation();

    StepStoppingCriterion stepStoppingCriterion_0 = 
      ((StepStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Steps"));

    stepStoppingCriterion_0.setMaximumNumberSteps(20);

    simulation_0.getSimulationIterator().run();

  }
}
