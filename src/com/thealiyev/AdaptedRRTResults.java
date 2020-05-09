package com.thealiyev;

public class AdaptedRRTResults {
    public static void main(String[] args) {
        GWOPathPlanningDynamicStationaryMethod gwoPathPlanningDynamicStationaryMethod = new GWOPathPlanningDynamicStationaryMethod();
        gwoPathPlanningDynamicStationaryMethod.GWO();

        IGWOPathPlanningDynamicStationaryMethod igwoPathPlanningDynamicStationaryMethod = new IGWOPathPlanningDynamicStationaryMethod();
        igwoPathPlanningDynamicStationaryMethod.GWO();

        ExGWOPathPlanningDynamicStationaryMethod exgwoPathPlanningDynamicStationaryMethod = new ExGWOPathPlanningDynamicStationaryMethod();
        exgwoPathPlanningDynamicStationaryMethod.GWO();
    }

}
