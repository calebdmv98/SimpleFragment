import java.io.File;
import java.util.Locale;

import org.hipparchus.ode.nonstiff.AdaptiveStepsizeIntegrator;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.util.FastMath;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;
import org.orekit.errors.OrekitException;
import org.orekit.forces.ForceModel;
import org.orekit.forces.drag.DragForce;
import org.orekit.forces.drag.DragSensitive;
import org.orekit.forces.drag.IsotropicDrag;
import org.orekit.forces.drag.atmosphere.Atmosphere;
import org.orekit.forces.drag.atmosphere.SimpleExponentialAtmosphere;
import org.orekit.forces.gravity.potential.GravityFieldFactory;
import org.orekit.forces.gravity.potential.NormalizedSphericalHarmonicsProvider;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.OrbitType;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.KeplerianPropagator;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.propagation.sampling.OrekitFixedStepHandler;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScale;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;

public class OrePractise {
	 
    public static void main(String[] args) {
        try {

        		//Set HDD a main Directory
            File home       = new File(System.getProperty("user.home"));
            //Search for OREkit folder in the main directory set above
            File orekitData = new File(home, "orekit-data");
            if (!orekitData.exists()) {
                System.err.format(Locale.US, "Failed to find %s folder%n",
                                  orekitData.getAbsolutePath());
                System.err.format(Locale.US, "You need to download %s from the %s page and unzip it in %s for this tutorial to work%n",
                                  "orekit-data.zip", "https://www.orekit.org/forge/projects/orekit/files",
                                  home.getAbsolutePath());
                System.exit(1);
            }
            
            DataProvidersManager manager = DataProvidersManager.getInstance();
            manager.addProvider(new DirectoryCrawler(orekitData));

            // Inertial frame 256
            Frame inertialFrame = FramesFactory.getEME2000();

            // Initial date in UTC time scale
            TimeScale utc = TimeScalesFactory.getUTC();
            AbsoluteDate initialDate = new AbsoluteDate(2017, 12, 20, 23, 30, 00.000, utc);

            // gravitation coefficient
            double mu =  3.986004415e+14;
            
            // Initial orbit parameters
            double a = 24396159; // semi major axis in meters
            double e = 0.72831215; // eccentricity
            double i = FastMath.toRadians(7); // inclination
            double omega = FastMath.toRadians(180); // perigee argument
            double raan = FastMath.toRadians(261); // right ascension of ascending node
            double lM = 0; // mean anomaly
            	
            // Orbit construction as Keplerian
            Orbit initialOrbit = new KeplerianOrbit(a, e, i, omega, raan, lM, PositionAngle.MEAN,
                                                    inertialFrame, initialDate, mu);
            //initial StatenDefinition
            SpacecraftState initialState = new SpacecraftState(initialOrbit);
            
            // Step integrator 
            final double minStep = 0.001;
            final double maxstep = 1000.0;
            final double positionTolerance = 10.0;
            final OrbitType propagationType = OrbitType.KEPLERIAN;
            final double[][] tolerances =
                    NumericalPropagator.tolerances(positionTolerance, initialOrbit, propagationType);
            AdaptiveStepsizeIntegrator integrator =
                    new DormandPrince853Integrator(minStep, maxstep, tolerances[0], tolerances[1]);

            // Propagator
            NumericalPropagator propagator = new NumericalPropagator(integrator);
            propagator.setOrbitType(propagationType);

            // Force Model
            final double ae = Constants.GRS80_EARTH_EQUATORIAL_RADIUS;
            OneAxisEllipsoid earth = new OneAxisEllipsoid(ae, Constants.WGS84_EARTH_FLATTENING, FramesFactory.getEME2000());
            earth.setAngularThreshold(1.e-6);
            Atmosphere atm = new SimpleExponentialAtmosphere(earth, 4.e-13, 500000.0, 60000.0);
            final double cd = 2.0;
            final double sf = 5.0;
            ForceModel dragNUM = new DragForce(atm, new IsotropicDrag(sf, cd));

            // Add force model to the propagator
            propagator.addForceModel(dragNUM);

            // Set up initial state in the propagator
            propagator.setInitialState(initialState);

            // Set up operating mode for the propagator as master mode
            // with fixed step and specialized step handler
            propagator.setMasterMode(1., new stepHandler());

            // Extrapolate from the initial to the final date
            SpacecraftState finalState = propagator.propagate(initialDate.shiftedBy(10.));
            KeplerianOrbit o = (KeplerianOrbit) OrbitType.KEPLERIAN.convertType(finalState.getOrbit());
            System.out.print("\nPVC Coords: " + o.initPVCoordinates()+"\n");         
            /*
            System.out.format(Locale.US, "Final state:%n%s %12.3f %10.8f %10.6f %10.6f %10.6f %10.6f%n",
                              finalState.getDate(),
                              o.getA(), o.getE(),
                              FastMath.toDegrees(o.getI()),
                              FastMath.toDegrees(o.getPerigeeArgument()),
                              FastMath.toDegrees(o.getRightAscensionOfAscendingNode()),
                              FastMath.toDegrees(o.getTrueAnomaly()));
     	
        } catch (OrekitException oe) {
            System.err.println(oe.getMessage());
        }
    }
    
    
    // Step Handler Class used to print on the output stream at the given step.
    private static class stepHandler implements OrekitFixedStepHandler {

        private stepHandler() {
            //private constructor
        }

        /*
        public  void init(final SpacecraftState s0, final AbsoluteDate t, final double step) {
            System.out.println("          date                a           e" +
                               "           i         \u03c9          \u03a9" +
                               "          \u03bd");
        }
        */
        public  void handleStep(SpacecraftState currentState, boolean isLast) {
            KeplerianOrbit o = (KeplerianOrbit) OrbitType.KEPLERIAN.convertType(currentState.getOrbit());
            /*
            System.out.format(Locale.US, "%s %12.3f %10.8f %10.6f %10.6f %10.6f %10.6f%n",
                              currentState.getDate(),
                              o.getA(), o.getE(),
                              FastMath.toDegrees(o.getI()),
                              FastMath.toDegrees(o.getPerigeeArgument()),
                              FastMath.toDegrees(o.getRightAscensionOfAscendingNode()),
                              FastMath.toDegrees(o.getTrueAnomaly()));
            */
            System.out.print("\nPVC Coords: " + o.initPVCoordinates()+"\n");      
            if (isLast) {
                System.out.println("this was the last step ");
                System.out.println();
            }
        }

    }

}

