/*
 * Copyright (c) 2023 ETH Zurich
 *
 * This file is part of remaster.
 *
 * remaster is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * remaster is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with remaster. If not, see <https://www.gnu.org/licenses/>.
 */

package remaster;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import remaster.reactionboxes.BDReactionBox;
import remaster.reactionboxes.ContinuousBDReactionBox;
import remaster.reactionboxes.PunctualBDReactionBox;

import java.io.PrintStream;
import java.util.*;

@Description("An object representing a stochastic birth-death trajectory." +
        "The birth-death model is specified via Functions representing" +
        "the various populations, and Reactions representing the reactions " +
        "producing the dynamics.  This object can be logged to produce a" +
        "TSV file which can be directly read into R for plotting.")
public class StochasticTrajectory extends AbstractBDTrajectory {

    public Input<Integer> maxRetriesInput = new Input<>("maxRetries",
            "Maximum number of times to retry simulation after failure " +
                    "to meet mustHave condition.",
            100);

    public Input<Double> logIntervalInput = new Input<>("logInterval",
            "Time interval between logging points. If <= 0 (default), logs every event.", -1.0);

    public Input<Double> tauLeapingIntervalInput = new Input<>("tauLeapingInterval",
            "Fixed time step for tau-leaping (days). If <= 0, uses exact Gillespie.", -1.0);

    public Input<Double> tauLeapingThresholdInput = new Input<>("tauLeapingThreshold",
            "Population threshold for hybrid tau-leaping. If > 0, uses exact SSA below threshold.", -1.0);

    public boolean currentTrajectoryValid;
    private boolean useTauLeaping = false;

    private double[] eventTimes;
    private double[] eventMultiplicities;
    private int[] eventReactionIdx;
    private int eventCount;
    private boolean eventsContainSamples;
    private Map<BDReactionBox, Integer> reactionBoxIndex;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        reactionBoxIndex = new HashMap<>();
        for (int i = 0; i < reactionBoxes.size(); i++)
            reactionBoxIndex.put(reactionBoxes.get(i), i);
        initEventStorage();

        double tau = tauLeapingIntervalInput.get();
        double threshold = tauLeapingThresholdInput.get();
        useTauLeaping = tau > 0.0;

        int retries = maxRetriesInput.get();
        while (retries >= 0 && !(currentTrajectoryValid =
                (useTauLeaping ? doTauLeapingSimulation(tau, threshold) : doSimulation()))) {
            retries -= 1;
            if (retries >= 0)
                Log.info("Trajectory simulation rejected: retrying");
            else
                Log.err.println("Failed to simulate trajectory satisfying " +
                        "mustHave condition. (maxRetries = " +
                        maxRetriesInput.get() + ")");
        }
    }

    @Override
    public boolean doSimulation() {
        state.resetToInitial();
        clearEvents();

        for (BDReactionBox reactionBox : reactionBoxes)
            reactionBox.resetInterval();

        List<BDReactionBox> reactionBoxesSortedByChangeTimes = new ArrayList<>(reactionBoxes);
        reactionBoxesSortedByChangeTimes.sort(Comparator.comparingDouble(BDReactionBox::getIntervalEndTime));

        double t = 0.0;

        while (true) {
            double a0 = 0.0;
            for (ContinuousBDReactionBox reactionBox : continuousReactionBoxes)
                a0 += reactionBox.updatePropensity();

            double delta = a0 == 0 ? Double.POSITIVE_INFINITY : Randomizer.nextExponential(a0);
            t += delta;

            BDReactionBox updatedReactionBox = reactionBoxesSortedByChangeTimes.get(0);
            if (maxTimeInput.get().getArrayValue() < updatedReactionBox.getIntervalEndTime()) {
                if (t > maxTimeInput.get().getArrayValue())
                    break;
            } else if (t > updatedReactionBox.getIntervalEndTime()) {
                t = updatedReactionBox.getIntervalEndTime();

                if (updatedReactionBox instanceof PunctualBDReactionBox) {
                    // Implement punctual reaction
                    double multiplicity = ((PunctualBDReactionBox) updatedReactionBox).implementEvent(true);
                    if (multiplicity < 0)
                        return false; // Occurs when n-events lack the necessary reactants to fire

                    if (multiplicity > 0)
                        addEvent(t, updatedReactionBox, multiplicity);
                }

                updatedReactionBox.incrementInterval();
                reactionBoxesSortedByChangeTimes
                        .sort(Comparator.comparingDouble(BDReactionBox::getIntervalEndTime));

                continue;
            }

            if (delta == Double.POSITIVE_INFINITY)
                break;

            double u = Randomizer.nextDouble() * a0;

            ContinuousBDReactionBox thisReactionBox = null;
            for (ContinuousBDReactionBox reaction : continuousReactionBoxes) {
                if (u < reaction.currentPropensity) {
                    thisReactionBox = reaction;
                    break;
                } else
                    u -= reaction.currentPropensity;
            }

            if (thisReactionBox == null)
                throw new IllegalStateException("Reaction selection loop fell through.");

            addEvent(t, thisReactionBox, 1);
            thisReactionBox.incrementState(state, 1);

            if (endCondition != null && endCondition.isMet()) {
                System.out.println("Trajectory termination condition met: " + endsWhenInput.get());
                break;
            }
        }

        if (acceptCondition != null && !acceptCondition.isMet()) {
            System.out.println("Trajectory acceptance condition not met: " + mustHaveInput.get());
            return false;
        }

        state.setFinal();
        return true;
    }

    private double getMaxVPopulation() {
        double maxV = 0.0;
        for (String popName : state.popIndices.keySet()) {
            if (popName.startsWith("V_")) {
                double val = state.get(popName, 0);
                if (val > maxV)
                    maxV = val;
            }
        }
        return maxV;
    }

    public boolean doTauLeapingSimulation(double tau, double threshold) {
        if (tau <= 0.0)
            return doSimulation();

        state.resetToInitial();
        clearEvents();

        for (BDReactionBox reactionBox : reactionBoxes)
            reactionBox.resetInterval();

        List<BDReactionBox> reactionBoxesSortedByChangeTimes = new ArrayList<>(reactionBoxes);
        reactionBoxesSortedByChangeTimes.sort(Comparator.comparingDouble(BDReactionBox::getIntervalEndTime));

        double t = 0.0;
        double maxTime = maxTimeInput.get().getArrayValue();

        while (true) {
            BDReactionBox updatedReactionBox = reactionBoxesSortedByChangeTimes.get(0);
            double nextChangeTime = updatedReactionBox.getIntervalEndTime();

            if (nextChangeTime <= t + 1e-12) {
                if (updatedReactionBox instanceof PunctualBDReactionBox) {
                    double multiplicity = ((PunctualBDReactionBox) updatedReactionBox).implementEvent(true);
                    if (multiplicity < 0)
                        return false;
                    if (multiplicity > 0)
                        addEvent(t, updatedReactionBox, multiplicity);
                }

                updatedReactionBox.incrementInterval();
                reactionBoxesSortedByChangeTimes
                        .sort(Comparator.comparingDouble(BDReactionBox::getIntervalEndTime));
                continue;
            }

            if (t >= maxTime)
                break;

            boolean useExact = threshold > 0.0 && getMaxVPopulation() < threshold;
            if (useExact) {
                double a0 = 0.0;
                for (ContinuousBDReactionBox reactionBox : continuousReactionBoxes)
                    a0 += reactionBox.updatePropensity();

                double delta = a0 == 0 ? Double.POSITIVE_INFINITY : Randomizer.nextExponential(a0);
                double nextStop = Math.min(nextChangeTime, maxTime);
                if (t + delta > nextStop) {
                    t = nextStop;
                    if (t >= maxTime)
                        break;
                    continue;
                }

                t += delta;
                double u = Randomizer.nextDouble() * a0;
                ContinuousBDReactionBox thisReactionBox = null;
                for (ContinuousBDReactionBox reaction : continuousReactionBoxes) {
                    if (u < reaction.currentPropensity) {
                        thisReactionBox = reaction;
                        break;
                    } else
                        u -= reaction.currentPropensity;
                }
                if (thisReactionBox == null)
                    throw new IllegalStateException("Reaction selection loop fell through.");

                addEvent(t, thisReactionBox, 1);
                thisReactionBox.incrementState(state, 1);
            } else {
                double migSum = 0.0;
                double[] migProps = new double[continuousReactionBoxes.size()];
                for (int i = 0; i < continuousReactionBoxes.size(); i++) {
                    ContinuousBDReactionBox reactionBox = continuousReactionBoxes.get(i);
                    if (reactionBox.isMigrationReaction()) {
                        double a = reactionBox.updatePropensity();
                        migProps[i] = a;
                        migSum += a;
                    }
                }

                double nextMigDelta = migSum == 0.0 ? Double.POSITIVE_INFINITY : Randomizer.nextExponential(migSum);

                double dt = tau;
                double nextStop = Math.min(nextChangeTime, maxTime);
                if (t + dt > nextStop)
                    dt = nextStop - t;
                if (t + nextMigDelta < t + dt)
                    dt = nextMigDelta;
                if (dt <= 0.0)
                    continue;

                double[] ks = new double[continuousReactionBoxes.size()];
                double[] props = new double[continuousReactionBoxes.size()];

                boolean anyPropensity = false;
                for (int i = 0; i < continuousReactionBoxes.size(); i++) {
                    ContinuousBDReactionBox reactionBox = continuousReactionBoxes.get(i);
                    if (reactionBox.isMigrationReaction())
                        continue;
                    double a = reactionBox.updatePropensity();
                    props[i] = a;
                    if (a > 0.0)
                        anyPropensity = true;
                }

                if (!anyPropensity && (Double.isInfinite(nextChangeTime) || nextChangeTime > maxTime) && nextMigDelta == Double.POSITIVE_INFINITY)
                    break;

                int retries = 0;
                final int maxRetries = 10;
                boolean validStep = false;
                double[] stateBackup = Arrays.copyOf(state.occupancies, state.occupancies.length);
                int eventCountBackup = eventCount;

                while (!validStep && retries <= maxRetries) {
                    for (int i = 0; i < continuousReactionBoxes.size(); i++) {
                        ContinuousBDReactionBox reactionBox = continuousReactionBoxes.get(i);
                        if (reactionBox.isMigrationReaction()) {
                            ks[i] = 0.0;
                            continue;
                        }
                        double lambda = props[i] * dt;
                        if (lambda <= 0.0) {
                            ks[i] = 0.0;
                            continue;
                        }

                        double k = Util.nextPoisson(lambda);
                        double maxCount = reactionBox.getMaxReactCount();
                        if (k > maxCount)
                            k = maxCount;
                        ks[i] = k;
                    }

                    // Apply and validate
                    for (int i = 0; i < continuousReactionBoxes.size(); i++) {
                        double k = ks[i];
                        if (k > 0.0) {
                            ContinuousBDReactionBox reactionBox = continuousReactionBoxes.get(i);
                            reactionBox.incrementState(state, k);
                        addEvent(t + dt, reactionBox, k);
                        }
                    }

                    if (state.isValid()) {
                        validStep = true;
                    } else {
                        // Roll back and retry with smaller step
                        System.arraycopy(stateBackup, 0, state.occupancies, 0, stateBackup.length);
                        eventCount = eventCountBackup;
                        dt *= 0.5;
                        retries += 1;
                    }
                }

                if (!validStep)
                    return false;

                t += dt;

                if (Math.abs(dt - nextMigDelta) < 1e-12 && migSum > 0.0) {
                    double u = Randomizer.nextDouble() * migSum;
                    ContinuousBDReactionBox migBox = null;
                    for (int i = 0; i < continuousReactionBoxes.size(); i++) {
                        if (migProps[i] <= 0.0)
                            continue;
                        if (u < migProps[i]) {
                            migBox = continuousReactionBoxes.get(i);
                            break;
                        } else
                            u -= migProps[i];
                    }
                    if (migBox != null) {
                        addEvent(t, migBox, 1);
                        migBox.incrementState(state, 1);
                    }
                }
            }

            if (endCondition != null && endCondition.isMet()) {
                System.out.println("Trajectory termination condition met: " + endsWhenInput.get());
                break;
            }
        }

        if (acceptCondition != null && !acceptCondition.isMet()) {
            System.out.println("Trajectory acceptance condition not met: " + mustHaveInput.get());
            return false;
        }

        state.setFinal();
        return true;
    }

    @Override
    public Node simulateTree() throws SimulationFailureException {
        if (!eventsContainSamples)
            throw new SimulationFailureException("Trajectory contains no samples.");

        if (!currentTrajectoryValid)
            throw new SimulationFailureException(
                    "Refusing to simulate tree from trajectory failing mustHave condition.");

        Map<ReactElement, List<Lineage>> lineages = new HashMap<>();

        state.resetToFinal();

        LineageFactory lineageFactory = new LineageFactory();

        for (int idx = eventCount - 1; idx >= 0; idx--) {
            BDReactionBox reactionBox = reactionBoxes.get(eventReactionIdx[idx]);
            double time = eventTimes[idx];
            long n = Math.round(eventMultiplicities[idx]);
            if (n <= 0)
                continue;

            if (!useTauLeaping || n == 1) {
                for (long i = 0; i < n; i++) {
                    reactionBox.incrementLineages(lineages, time,
                            lineageFactory, false);
                    reactionBox.incrementState(state, -1);
                }
                continue;
            }

            double pInclude = reactionBox.getLineageInclusionProbability(lineages);
            long k;
            if (pInclude <= 0.0) {
                k = 0;
            } else if (pInclude >= 1.0) {
                k = n;
            } else {
                k = Math.round(Util.nextBinomial(n, pInclude));
            }
            reactionBox.incrementState(state, -n);
            for (long i = 0; i < k; i++) {
                reactionBox.incrementLineages(lineages, time,
                        lineageFactory, true);
            }
        }

        List<Lineage> rootLineages = new ArrayList<>();
        for (ReactElement pop : lineages.keySet())
            rootLineages.addAll(lineages.get(pop));

        // Restrict reactions to ensure this is impossible
        if (rootLineages.isEmpty())
            throw new SimulationFailureException("No lineages remaining.");

        // Might allow this in future
        if (rootLineages.size() > 1)
            throw new SimulationFailureException("Multiple lineages remaining.");

        lineageFactory.numberInternals(rootLineages.get(0));
        lineageFactory.computeAgesFromTimes(rootLineages.get(0));

        return rootLineages.get(0);
    }

    @Override
    public void log(long sample, PrintStream out) {
        state.resetToInitial();

        double logInterval = logIntervalInput.get();

        if (logInterval <= 0) {
            // Original behavior: log every event
            state.addToLog(out, sample, 0, true);

            for (int i = 0; i < eventCount; i++) {
                BDReactionBox reactionBox = reactionBoxes.get(eventReactionIdx[i]);
                reactionBox.incrementState(state, eventMultiplicities[i]);
                state.addToLog(out, sample, eventTimes[i], false);
            }
        } else {
            // Downsampled logging at fixed intervals
            state.addToLog(out, sample, 0, true);
            double nextLogTime = logInterval;

            for (int i = 0; i < eventCount; i++) {
                double eventTime = eventTimes[i];
                BDReactionBox reactionBox = reactionBoxes.get(eventReactionIdx[i]);
                while (nextLogTime < eventTime) {
                    // Log the *current* state (valid from previous event time up to event.time)
                    state.addToLog(out, sample, nextLogTime, false);
                    nextLogTime += logInterval;
                }
                reactionBox.incrementState(state, eventMultiplicities[i]);
            }

            // Log any remaining intervals up to the last event time (or end of simulation)
            // Note: The loop naturally handles everything *before* the last event.
            // The state after the last event corresponds to the final state.
            // If nextLogTime is still within simulation bounds (e.g. if we consider
            // maxTime), we might log more.
            // For now, we stop at the last event to remain consistent with original
            // behavior which stops at last event time.
            // Typically, one might want to include the final state at t=maxTime if needed.
            if (eventCount > 0 && nextLogTime <= eventTimes[eventCount - 1]) {
                state.addToLog(out, sample, nextLogTime, false);
            }
        }

        out.print("\t");
    }

    private void initEventStorage() {
        eventTimes = new double[1024];
        eventMultiplicities = new double[1024];
        eventReactionIdx = new int[1024];
        eventCount = 0;
        eventsContainSamples = false;
    }

    private void clearEvents() {
        eventCount = 0;
        eventsContainSamples = false;
    }

    private void addEvent(double time, BDReactionBox reactionBox, double multiplicity) {
        if (multiplicity <= 0.0)
            return;
        forceAddEvent(time, reactionBox, multiplicity);
    }

    private void forceAddEvent(double time, BDReactionBox reactionBox, double multiplicity) {
        if (multiplicity <= 0.0)
            return;
        if (eventCount == Integer.MAX_VALUE) {
            throw new IllegalStateException("Event count exceeded integer limit. "
                    + "Simulation generated too many events; increase scale_factor or tau_leaping_dt.");
        }
        if (eventCount >= eventTimes.length) {
            if (eventTimes.length > Integer.MAX_VALUE / 2) {
                throw new IllegalStateException("Event storage exceeded Java array limits. "
                        + "Simulation generated too many events; increase scale_factor or tau_leaping_dt.");
            }
            int newSize = eventTimes.length * 2;
            eventTimes = Arrays.copyOf(eventTimes, newSize);
            eventMultiplicities = Arrays.copyOf(eventMultiplicities, newSize);
            eventReactionIdx = Arrays.copyOf(eventReactionIdx, newSize);
        }
        Integer idx = reactionBoxIndex.get(reactionBox);
        if (idx == null)
            throw new IllegalStateException("Reaction box index not found for event.");
        eventTimes[eventCount] = time;
        eventMultiplicities[eventCount] = multiplicity;
        eventReactionIdx[eventCount] = idx;
        eventCount += 1;
        if (reactionBox.producesSamples)
            eventsContainSamples = true;
    }

    public int getEventCount() {
        return eventCount;
    }

}
