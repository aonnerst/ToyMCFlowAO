const int NC = 3; // # of centralities
TString strCentrality[NC] = {"0-5\%","20-30\%","50-60\%"};
const float inputNch[NC] ={1943, 786, 183};
const int NH = 6; // # of harmonic orders
//rows: input vn from n=2 to n=7
//cols: vn for multiplicities from central to peripheral
const float inputVn[NH][NC] = {
    {0.044325430932513066,0.13409520177544662,0.11868377868652139},
    {0.020319043993158488,0.06147004658105161,0.05440535759432017},
    {0.009314371910439144,0.028178238867966162,0.02493974296844666},
    {0.0042697640752778244,0.01291707408507019,0.011432528097143131},
    {0.001957285518962443,0.005921264409071094,0.0052407396041462785},
    {0.0008972314477330478,0.0027143431996458607,0.0024023865382260113},
};