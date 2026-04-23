#ifndef COUNTY_H_
#define COUNTY_H_

#include <cstdint>

using namespace std;

class County
{
public:
	enum CountyArchetypes {
		
	};

	// TEMPLATE? or just copy it for each one?
	struct CountyArchetypePreset {

	};

	struct CountyPopulation {
		uint64_t populationSize;
		// 0 - 1
		float populationDensity;

	};
	
	struct CountyLimits {
		// How many patients at the same time this county can handle
		uint64_t startingHealthcareCapacity;
		
		// How many extra patients can be temporarily taken care of
		// beyond capacity
		uint64_t reserveCapacity;

		// How long the extra patients can be taken care of
		int reverseLastingTimesteps;

	};

	//1111111111111111111111111111111111111111111111111111111111111111111111111111
	//
		// How well the county can absorb disruption before instability starts rising rapidly
		economicBase
		// How well a county maintains its resources
		// in case of scarcity and logistical issues
		supplyPreservationAbility
		// How connected it is to its neighbors
		// Index of how much tourists visit from neighboring counties
		localTravelConnection
		// How connected it is to non-neighboring counties
		// Index of how much tourists visit from non-neighboring counties
		globalTravelConnection
		// Baseline trust of the population in local institutions
		const baselineTrust
		
	//2222222222222222222222222222222222222222222222222222222222222222222222222222
	// TRUE VALUES
		// Actual number of infected people
		trueInfectedPopulation
		// Current spread potential
		contagiousPressure
		// How much of the current healthcare system is currently being used up
		hospitalLoad
		// How many hospitalizations the county can actually support (can change)
		// Based off starting healthcare capacity, changes based on other factors
		effectiveHospitalCapacity
		// Index of how many public interactions still happen, such as gatherings
		publicInteractionLevel
		// Index of how much people obey the policies set in place
		complianceLevel
		// How worn down the population is from prolonged restrictions, 
		// crisis, or repeated emergency measures put in place
		// Increases risk of backlash, negative events
		// shifts balance between slight negative events and disastrous ones
		fatigueLevel
		// How much the population of the county trusts local institutions
		// Based off starting trust level, changes based on other factors
		publicTrust
		// How socially unstable the county currently is
		// Protests, disorder, civil pushback, breakdown, rioting
		// Higher levels make control more difficult and interferes
		// with healthcare, reporting to the coordinator and policy enforcement
		unrestLevel
		// How functional the county’s economy and daily-life support system currently are
		// Low economic stability can cause disruptions to be longer, 
		// make it harder to sustain the restrictions, increase mortality due to starvation
		economicStability
		// How intact the supply chain of this county is
		// Worse condition means less medical supplies, increased fatigue,
		// less food, worse economic stability
		supplyCondition
		// How accurate and recent the reports are for this county
		// Bad reporting quality = miscounted cases, understated burden,
		// stale information
		localReportingQuality
		// How accurate the received global reports are
		// Poor communication may cause staff to misreport information
		// received from the global reporting
		// If it is really poor, global reporting may be missed entirely 
		globalReceivedReportQuality
		// Accumulated burden of the recent policy actions
		// Lingering weight of prior interventions
		// Decreases over time
		// Fatigue, trust, and unrest are hurt by high policy burden
		policyBurden
		// How much policy burden the public forgives/forgets every timestep
		policyBurdenDecreasePerTimestep
		// How much death is impacting the population
		// more death -> worse public sanity factors
		// even if the infection count is low
		deathBurden
		cumulitiveDeaths
		recentDeathRate

	//3333333333333333333333333333333333333333333333333333333333333333333333333333
		// Officially visible infection count
		reportedCases
		// Officially visible hospital burden
		reportedHospitalBurden
		// Officially visible public interaction level
		estimatedPublicInteractionLevel
		// Enum?
		// Officially visible indicator of how bad the social instability is
		// Reported as a signal based on events and true data:
		// Protests, Public Disorder, Disobedience, Tensions, Riots, Backlash,
		// System Collapse
		unrestSignal
		// Struct?
		// Officially visible shortages, logistical problems, resource issues
		supplyReport
		// Estimate of how much the public trusts the decision-maker
		trustEstimate
		// Struct?
		// Reported cases rising, hospital burden changing, unrest changing
		// Derived variable
		trendEstimate
		// Struct?
		// Last global report available
		// Can indicate national alert level, global events / issues
		globalReport
		// Struct?
		// Similar to global report but is much more frequent
		// Sent by neighbors to coordinate 
		neighborReport
		
	//4444444444444444444444444444444444444444444444444444444444444444444444444444
		// How much ordinary people are pushing the policy in their preferred direction
		// Rises when social instability rises
		generalPopulationPressure
		// How much vulnerable people are pushing the policy in their preferred direction
		vulnerablePopulationPressure
		// How much essential workers are pushing the policy in their preferred direction
		essentialWorkerPopulationPressure
		// How much hospital owners are pushing the policy in their preferred direction
		hospitalOwnerPopulationPressure
		// How much politicians are pushing the policy in their preferred direction
		politicianPopulationPressure
		// How much economic stakeholders are pushing the policy in their preferred direction
		economicStakeholderPopulationPressure

private:

};

#endif


/*

Paper Worthy
Internal AI competition inside each county, where multiple strategic actors actively fight for policy control instead of only contributing pressure weights
Coalition and takeover logic that determines which internal actor or bloc controls county policy at a given timestep
Scenario classes with genuinely different optimal strategies, so that one controller style is not dominant across all runs
Robust comparison framework between heuristic controllers, planner-based controllers, and learned controllers under the same seeded scenarios
County identity evolution over time, where counties can become resilient, brittle, radicalized, exhausted, stabilized, hollowed-out, or politically captured
Systemic breakdown cascades where one county’s collapse can trigger multi-county failure chains in healthcare, logistics, reporting, or governance
Hidden regime shifts that force controllers to infer changes in world behavior, such as variant changes, reporting regime changes, social fatigue shifts, or logistics regime shifts
Events that directly reshape county control, legitimacy, or ruling coalition rather than only modifying scalar state variables
Cross-county strategic or ideological blocs that influence policy alignment and coordinated reaction across multiple counties
Long-memory political and social state beyond simple policy burden, such as historical betrayal, chronic instability, hardening of public resistance, or accumulated institutional damage

Stretch
Dynamic controller profile switching inside a county without full internal AI competition
Long event cascades with persistent event memory and delayed aftereffects
Region-to-region resource transfer systems such as support sharing, reserve movement, or emergency aid
Explicit legitimacy and governance collapse mechanics
Long-lasting policy shadow effects that reshape future county behavior even after the immediate policy is gone
Long-range travel network overlay beyond nearest-neighbor spread
Contradictory, manipulated, degraded, or suppressed reporting
County archetype preset system with structured starting templates
Persistent global regime state that changes the behavior of all counties at once
Inter-county unrest or political contagion beyond direct disease spread

*/
