//
// Cattani-et-al-C++code
//
//  Created by Anna Cattani on 4/13/23.
//
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream> //this is the library needed for file streaming
#include <stdlib.h>
#include <time.h> // needed to seed srand
#include <cmath>
#include <string>
#include <unistd.h>
#include <vector>

/*-- constants --*/

#define  SVCOUNT     67 // number of state variables
#define  SVCOUNT2    10 // number of each state variables

#define  VIPCOUNT   3  // number of VIP cells
#define  SOMCOUNT   3  // number of SOM cells
#define  PVCOUNT    3  // number of PV cells
#define  ECSCOUNT   10 // number od ECS cells
#define  PNCOUNT    10 // number of PN cells
#define  PYRCSCOUNT 1  // number of auxiliary CS excitatory neurons
#define  PYRUSCOUNT 1  // number of auxiliary US excitatory neurons
#define  CURRCOUNT  3

#define  START        0
#define  END        40000 //ms
#define  TIMESTEP   0.05

#define trialCOUNT 1


/*-- parameters that should be set for each simulation --*/
int ext_stimul = 3; // 0: neither CS nor US; 1: only CS; 2: only US; 3: US+CS
int condit = 1; // 0: before FC; 1: during FC; 2: after FC;
int exper = 0; //0: PV and SOM; 1: no SOM; 2: no PV; 3: neither PV nor SOM; 4: only PV+PN connections; 5: no connections at all; 6: no VIP; 7: PV-only; 8: SOM-only
int stdp_rule = 1; //1: asym with larger tau and same amplitude;


/*-- connection strengths --*/
double gGABA_fromVIPtoSOM = 1.0/VIPCOUNT;
double gGABA_fromVIPtoPV = 1.0/VIPCOUNT;
double gGABA_fromSOMtoECS = 0.4/SOMCOUNT;
double gGABA_fromSOMtoPN = 0.4/SOMCOUNT;
double gGABA_fromPVtoSOM = 0;
double gGABA_fromPVtoPN = 0.5/PVCOUNT;
double gAMPA_fromPNtoPV = 0.5;
double gGABA_fromPVtoECS = 0.4/PVCOUNT;
double gAMPA_fromPNtoVIP[VIPCOUNT][PNCOUNT];
double gAMPA_fromPYRCStoECS[ECSCOUNT];
double gAMPA_fromPYRUStoPN[PNCOUNT];
double gAMPA_fromPYRCStoPV;

double gAMPA_fromECStoPN[ECSCOUNT];


/*-- plasticity parameters --*/
double tau_plus;
double tau_minus;
double A_plus;
double A_minus;
double tau_plus_VIP;
double tau_minus_VIP;
double A_plus_VIP;
double A_minus_VIP;

/*-- external inputs parameters --*/
double IappVIP[VIPCOUNT];
double IappPN = 0.35;
double IappECS = 0.45;
double strength_poisson_CS = 30;
double strength_poisson_PN = 30;


/*-- function prototypes --*/
float box_muller(float m, float s);    /* normal random variate generator */
int IntRK4( int svCount, int svCount2, double sv[SVCOUNT][SVCOUNT2], double times, double dt, double LFP[0], double I[CURRCOUNT], double gAMPA_fromECStoPN[ECSCOUNT], double gAMPA_fromPYRCStoPV, double isPoissonSpikeCS, double isPoissonSpikePN);
void Derive( int svCount, int svCount2, double y[SVCOUNT][SVCOUNT2], double k[SVCOUNT][SVCOUNT2], double times, double LFP[0],double I[CURRCOUNT], double gAMPA_fromECStoPN[ECSCOUNT], double gAMPA_fromPYRCStoPV, double isPoissonSpikeCS, double isPoissonSpikePN);
double rnd();

/*-- main program --*/
using namespace std;
int num_steps = double(END)/TIMESTEP;

int main( void )
{
    if (exper == 1) // no SOM
    {
        gGABA_fromSOMtoECS = 0;
        gGABA_fromSOMtoPN = 0;
    }
    if (exper == 2) // no PV
    {
        gGABA_fromPVtoPN = 0;
        gAMPA_fromPNtoPV = 0;
        gGABA_fromPVtoECS = 0;
    }
    if (exper == 3) // VIP-only
    {
        gGABA_fromSOMtoECS = 0;
        gGABA_fromSOMtoPN = 0;
        gGABA_fromPVtoPN = 0;
        gAMPA_fromPNtoPV = 0;
        gGABA_fromPVtoECS = 0;
    }
    if (exper == 6) // no VIP
    {
        gGABA_fromVIPtoSOM = 0;
        gGABA_fromVIPtoPV = 0;
    }
    if (exper == 7) // PV-only
    {
        gGABA_fromVIPtoSOM = 0;
        gGABA_fromVIPtoPV = 0;
        gGABA_fromSOMtoECS = 0;
        gGABA_fromSOMtoPN = 0;
    }

    gAMPA_fromPYRCStoPV = 0.2;
    gAMPA_fromPYRCStoECS[0] = 0.2;
    int z;
    for (z=1; z<ECSCOUNT; z++) {
        gAMPA_fromPYRCStoECS[z] = 0; //only ECS[0] receives CS
    }
    gAMPA_fromPYRUStoPN[0] = 0.2;
    for (z=1; z<PNCOUNT; z++) {
        gAMPA_fromPYRUStoPN[z] = 0; //only PN[0] receives US
    }

    //ext_stimul = 0: neither CS nor US; 1: only CS; 2: only US; 3: US+CS
    if (ext_stimul == 0) //0: neither CS nor US;
    {
        IappVIP[0] = 4.5;
        IappVIP[1] = 4;
        IappVIP[2] = 3.5;

        strength_poisson_PN = 0;
        strength_poisson_CS = 0;
    }

    if (ext_stimul == 1) // only CS;
    {
        IappVIP[0] = 4.1;
        IappVIP[1] = 4;
        IappVIP[2] = 3.9;

        strength_poisson_PN = 0;
    }

    if (ext_stimul == 2) // only US;
    {
        IappVIP[0] = 5;
        IappVIP[1] = 5;
        IappVIP[2] = 5;
        IappPN = 0.5;

        strength_poisson_CS = 0;
    }

    if (ext_stimul == 3) // 3: US+CS
    {
        IappVIP[0] = 5;
        IappVIP[1] = 5;
        IappVIP[2] = 5;
        IappPN = 0.5;
    }

    tau_plus = 14;
    tau_minus = 28;
    A_plus = 5e-3;
    A_minus = 5e-3;

    tau_plus_VIP = 14;
    tau_minus_VIP = 28;
    A_plus_VIP = 0.65e-3;
    A_minus_VIP = 0.3e-3;

    //******//

    int trial;
    for (trial=0; trial<trialCOUNT; trial++)
    {
        if (condit == 2) // 2: after FC
        {
            int z;
            gAMPA_fromECStoPN[0] = 0.18;
            for (z=1; z<ECSCOUNT; z++) {
                gAMPA_fromECStoPN[z] = 0.0001;
            }
            gAMPA_fromPNtoVIP[0][0] = 0.04;
            gAMPA_fromPNtoVIP[1][0] = 0.04;
            gAMPA_fromPNtoVIP[2][0] = 0.04;
            for (z=1; z<PNCOUNT; z++){
                gAMPA_fromPNtoVIP[0][z] = 0.01;
                gAMPA_fromPNtoVIP[1][z] = 0.01;
                gAMPA_fromPNtoVIP[2][z] = 0.01;
            }
        }
        else // before FC
        {
            int z;
            gAMPA_fromECStoPN[0] = 0.0001;
            for (z=0; z<ECSCOUNT; z++)
            {
                gAMPA_fromECStoPN[z] = 0.0001;
            }
            for (z=0; z<PNCOUNT; z++){
                gAMPA_fromPNtoVIP[0][z] = 0.01;
                gAMPA_fromPNtoVIP[1][z] = 0.01;
                gAMPA_fromPNtoVIP[2][z] = 0.01;
            }
        }

        /* Output to file */

        std::string experiment = std::to_string(exper);
        std::string external_stimulation = std::to_string(ext_stimul);
        std::string condition = std::to_string(condit);
        std::string STDP_rule = std::to_string(stdp_rule);
        std::string tmp_trial = std::to_string(trial);

        string path_("/Users/acattani/Documents/data/fearconditioning/heterogeneous/check/");

        string s1("STDP_");
        string results = path_ + "cond" +  condition + "_stim" + external_stimulation + "_exp" + experiment + "_rule" + STDP_rule + s1 + tmp_trial + ".txt";
        ofstream outputFile(results, ios::out);

        string s2("param_");
        string results2 = path_ + "cond" +  condition + "_stim" + external_stimulation + "_exp" + experiment + "_rule" + STDP_rule + s2 + tmp_trial + ".txt";
        ofstream outputFile2(results2, ios::out);

        string s6("AMPA_ECS_PN_");
        string results6 = path_ + "cond" +  condition + "_stim" + external_stimulation + "_exp" + experiment + "_rule" + STDP_rule + s6 + tmp_trial + ".txt";
        ofstream outputFile6(results6, ios::out);

        string s7("LFP_");
        string results7 = path_ + "cond" +  condition + "_stim" + external_stimulation + "_exp" + experiment + "_rule" + STDP_rule + s7 + tmp_trial + ".txt";
        ofstream outputFile7(results7, ios::out);



        if (!outputFile)
        {
            cerr << "File could not be opened" << endl;
            exit(1);
        }

        /* Seed random number generator */
        srand ((unsigned int) time(NULL));


        /*** define variables ****/
        double sv[ SVCOUNT ][SVCOUNT2];
        double LFP[0];
        double I[3];
        double PN_ispk[PNCOUNT]; double PN_ispktime[PNCOUNT]; double PN_svp[PNCOUNT];
        double ECS_ispk[ECSCOUNT]; double ECS_ispktime[ECSCOUNT]; double ECS_svp[ECSCOUNT];
        double PV_ispk; double PV_ispktime; double PV_svp[PVCOUNT];
        double PYRCS_ispk; double PYRCS_ispktime; double PYRCS_svp;
        double VIP_ispk[VIPCOUNT]; double VIP_svp[VIPCOUNT]; double VIP_ispktime[VIPCOUNT];

        double times; int itime; int itime2; int ntime;

        int z;

        /*** Initial conditions for VIP cells ****/
        for (z=0; z<VIPCOUNT; z++)
        {
            sv[0][z] = ((-64 + 66) * ((float)rand() / RAND_MAX)) - 66;
            sv[1][z] = 0.32*(sv[0][z]+54)/(1-exp(-(sv[0][z]+54)/4)) / (0.32*(sv[0][z]+54)/(1-exp(-(sv[0][z]+54)/4))  + 0.28*(sv[0][z]+27)/(exp((sv[0][z]+27)/5)-1));
            sv[2][z] = 0.128*exp(-(sv[0][z]+50)/18) / (0.128*exp(-(sv[0][z]+50)/18)  +  4/(1+exp(-(sv[0][z]+27)/5)));
            sv[3][z] = 0.032*(sv[0][z]+52)/(1-exp(-(sv[0][z]+52)/5)) / (0.032*(sv[0][z]+52)/(1-exp(-(sv[0][z]+52)/5))   +  0.5*exp(-(sv[0][z]+57)/40));
            sv[5][z] = 1/(1 + exp(-(sv[0][z]+50)/20));
            sv[6][z] = 1/(1 + exp((sv[0][z]+70)/6));

            sv[7][z] = 0.0001;
            sv[8][z] = 0.0001;

            VIP_svp[z] = 0; VIP_ispk[z] = 0; VIP_ispktime[z] = 0;
        }

        /*** Initial conditions for SOM cells ****/
        for (z=0; z<SOMCOUNT; z++)
        {
            sv[10][z] = ((-60 + 65) * ((float)rand() / RAND_MAX)) - 65;
            sv[11][z] = (-0.1 * (sv[10][z] + 23) / (exp(-0.1 * (sv[10][z] + 23)) - 1)) / (-0.1 * (sv[10][z] + 23) / (exp(-0.1 * (sv[10][z] + 23)) - 1) + 4 * exp(-(sv[10][z] + 48) / 18));
            sv[12][z] = (0.07 * exp(-(sv[10][z] + 37) / 20)) / (0.07 * exp(-(sv[10][z] + 37) / 20) + 1 / (exp(-0.1 * (sv[10][z] + 7)) + 1));
            sv[13][z] = (-0.01 * (sv[10][z] + 27) / (exp(-0.1 * (sv[10][z] + 27)) - 1)) / (-0.01 * (sv[10][z] + 27) / (exp(-0.1 * (sv[10][z] + 27)) - 1) + 0.125 * exp(-(sv[10][z] + 37) / 80));
            sv[14][z] = 1 / (1 + exp((sv[10][z] + 79.2) / 9.78)); //H current: (h_f_inf - y_O[4]) / tauh_f
            sv[15][z] = pow(1 / (1 + exp((sv[10][z] + 2.83) / 15.9)), 58);
            sv[16][z] =  1 / (1 + exp(-(sv[10][z] + 38) / 6.5)); // P current (persistent sodium current): (pinf - y_O[6]) / taup

            sv[17][z] = 0.0001; // GABA gating variable
        }

        /**** Initial conditions for PV cells ****/
        for (z=0; z<PVCOUNT; z++)
        {
            sv[20][z] = ((-60 + 65) * ((float)rand() / RAND_MAX)) - 65;
            sv[21][z] = 0.32*(sv[20][z]+54)/(1-exp(-(sv[20][z]+54)/4)) / (0.32*(sv[20][z]+54)/(1-exp(-(sv[20][z]+54)/4))  + 0.28*(sv[20][z]+27)/(exp((sv[20][z]+27)/5)-1));
            sv[22][z] = 0.128*exp(-(sv[20][z]+50)/18) / (0.128*exp(-(sv[20][z]+50)/18)  +  4/(1+exp(-(sv[20][z]+27)/5)));
            sv[23][z] = 0.032*(sv[20][z]+52)/(1-exp(-(sv[20][z]+52)/5)) / (0.032*(sv[20][z]+52)/(1-exp(-(sv[20][z]+52)/5))   +  0.5*exp(-(sv[20][z]+57)/40));
            sv[24][z] = 0.0001;

            PV_svp[z] = 0;
        }


        /**** Initial conditions for ECS cells ****/
        for (z=0; z<ECSCOUNT; z++)
        {
            sv[30][z] = ((-60 + 65) * ((float)rand() / RAND_MAX)) - 65;
            sv[31][z] = (0.1 * (sv[30][z] + 35) / (1 - exp(-(sv[30][z] + 35) / 10))) / (0.1 * (sv[30][z] + 35) / (1 - exp(-(sv[30][z] + 35) / 10)) + 4 * exp(-(sv[30][z] + 60) / 18));
            sv[32][z] = (0.07 * exp(-(sv[30][z] + 58) / 20)) / (0.07 * exp(-(sv[30][z] + 58) / 20) + 1. / (1 + exp(-(sv[30][z] + 28) / 10)));
            sv[33][z] = (-0.01 * (sv[30][z] + 34) / (-1 + exp(-(sv[30][z] + 34) / 10))) / (-0.01 * (sv[30][z] + 34) / (-1 + exp(-(sv[30][z] + 34) / 10)) +  0.125 * exp(-(sv[30][z] + 44) / 80));
            sv[34][z] =1. / (1 + exp(-(sv[30][z] + 35) / 10));
            sv[35][z] = 0.0001;

            ECS_svp[z] = 0;
        }

        /**** Initial conditions for PN cells ****/
        for (z=0; z<PNCOUNT; z++)
        {
            sv[40][z] = ((-60 + 65) * ((float)rand() / RAND_MAX)) - 65;
            sv[41][z] = (0.1 * (sv[40][z] + 35) / (1 - exp(-(sv[40][z] + 35) / 10))) / (0.1 * (sv[40][z] + 35) / (1 - exp(-(sv[40][z] + 35) / 10)) + 4 * exp(-(sv[40][z] + 60) / 18));
            sv[42][z] = (0.07 * exp(-(sv[40][z] + 58) / 20)) / (0.07 * exp(-(sv[40][z] + 58) / 20) + 1. / (1 + exp(-(sv[40][z] + 28) / 10)));
            sv[43][z] = (-0.01 * (sv[40][z] + 34) / (-1 + exp(-(sv[40][z] + 34) / 10))) / (-0.01 * (sv[40][z] + 34) / (-1 + exp(-(sv[40][z] + 34) / 10)) +  0.125 * exp(-(sv[40][z] + 44) / 80));
            sv[44][z] =1. / (1 + exp(-(sv[40][z] + 35) / 10));
            sv[45][z] = 0.0001;

            PN_svp[z] = 0;

            sv[46][z] = 0; //PN_Mfun
            sv[47][z] = 0; //PN_Pfun (this is for plasticity F -> VIP
        }

        /**** Initial conditions for PYR CS cell ****/
        for (z=0; z<PYRCSCOUNT; z++)
        {
            sv[50][z] = ((-60 + 65) * ((float)rand() / RAND_MAX)) - 65;
            sv[51][z] = (0.1 * (sv[50][z] + 35) / (1 - exp(-(sv[50][z] + 35) / 10))) / (0.1 * (sv[50][z] + 35) / (1 - exp(-(sv[50][z] + 35) / 10)) + 4 * exp(-(sv[50][z] + 60) / 18));
            sv[52][z] = (0.07 * exp(-(sv[50][z] + 58) / 20)) / (0.07 * exp(-(sv[50][z] + 58) / 20) + 1. / (1 + exp(-(sv[50][z] + 28) / 10)));
            sv[53][z] = (-0.01 * (sv[50][z] + 34) / (-1 + exp(-(sv[50][z] + 34) / 10))) / (-0.01 * (sv[50][z] + 34) / (-1 + exp(-(sv[50][z] + 34) / 10)) +  0.125 * exp(-(sv[50][z] + 44) / 80));
            sv[54][z] = 1. / (1 + exp(-(sv[50][z] + 35) / 10));
            sv[55][z] = 0.0001;
        }

        for (z=0; z<PYRUSCOUNT; z++)
        {
            sv[60][z] = ((-60 + 65) * ((float)rand() / RAND_MAX)) - 65;
            sv[61][z] = (0.1 * (sv[60][z] + 35) / (1 - exp(-(sv[60][z] + 35) / 10))) / (0.1 * (sv[60][z] + 35) / (1 - exp(-(sv[60][z] + 35) / 10)) + 4 * exp(-(sv[60][z] + 60) / 18));
            sv[62][z] = (0.07 * exp(-(sv[60][z] + 58) / 20)) / (0.07 * exp(-(sv[60][z] + 58) / 20) + 1. / (1 + exp(-(sv[60][z] + 28) / 10)));
            sv[63][z] = (-0.01 * (sv[60][z] + 34) / (-1 + exp(-(sv[60][z] + 34) / 10))) / (-0.01 * (sv[60][z] + 34) / (-1 + exp(-(sv[60][z] + 34) / 10)) +  0.125 * exp(-(sv[60][z] + 44) / 80));
            sv[64][z] = 1. / (1 + exp(-(sv[60][z] + 35) / 10));
            sv[65][z] = 0.0001;
        }


        for (z=0; z<PNCOUNT; z++) {
            PN_ispktime[z] = -1;
        }
        for (z=0; z<ECSCOUNT; z++) {
            ECS_ispktime[z] = -1;
        }
        for (z=0; z<VIPCOUNT; z++) {
            VIP_ispktime[z] = -1;
        }

        PV_ispktime = -1;
        PYRCS_ispktime = -1;


        double isPoissonSpikeCS;
        double isPoissonSpikePN;
        isPoissonSpikeCS = 0;
        isPoissonSpikePN = 0;

        double lambda = 800; //per second
        double p = lambda * TIMESTEP/1000;


        outputFile2<<gGABA_fromVIPtoSOM<<","<<gGABA_fromVIPtoPV<<","<<gGABA_fromSOMtoECS<<","<<gGABA_fromSOMtoPN<<","<<gGABA_fromPVtoPN<<","<<gAMPA_fromPNtoPV<<","<<gGABA_fromPVtoECS<<","<<gAMPA_fromPNtoVIP[0][0]<<","<<gAMPA_fromPNtoVIP[1][0]<<","<<gAMPA_fromPNtoVIP[2][0]<<","<<gAMPA_fromPNtoVIP[0][1]<<","<<gAMPA_fromPYRCStoECS[0]<<","<<gAMPA_fromPYRCStoPV<<","<<gAMPA_fromPYRUStoPN<<endl; // save the parameters used

        times = START;
        /* start main time loop - loop until stopping time reached */
        while ( times <  END ) {
            ntime = times; // makes time an integer
            itime = ntime % 5;  // itime is an integer mod 5 --> (fs = 200), mod2 --> fs = 500

            IntRK4(SVCOUNT, SVCOUNT2, sv, times, TIMESTEP, LFP, I, gAMPA_fromECStoPN, gAMPA_fromPYRCStoPV,
                   isPoissonSpikeCS, isPoissonSpikePN);

            // generate the CS (isPoissonSpikeCS) and US (isPoissonSpikePN)
            isPoissonSpikeCS = 0;
            isPoissonSpikePN = 0;
            double randomnumber = rand() / double(RAND_MAX);
            if (randomnumber < p) {
                isPoissonSpikeCS = 1;
            }
            double randomnumber2 = rand() / double(RAND_MAX);
            if (randomnumber2 < p) {
                isPoissonSpikePN = 1;
            }

            //* file to be saved *//
            outputFile<<times<<","<<sv[0][0]<<","<<sv[0][1]<<","<<sv[0][2]<<","<<sv[10][0]<<","<<sv[10][1]<<","<<sv[10][2]<<","<<sv[20][0]<<","<<sv[30][0]<<","<<sv[30][1]<<","<<sv[30][2]<<","<<sv[30][3]<<","<<sv[30][4]<<","<<sv[30][5]<<","<<sv[30][6]<<","<<sv[30][7]<<","<<sv[30][8]<<","<<sv[30][9]<<","<<sv[40][0]<<","<<sv[40][1]<<","<<sv[40][2]<<","<<sv[40][3]<<","<<sv[40][4]<<","<<sv[40][5]<<","<<sv[40][6]<<","<<sv[40][7]<<","<<sv[40][8]<<","<<sv[40][9]<<","<<sv[50][0]<<","<<sv[60][0]<<","<<sv[36][0]<<","<<sv[46][0]<<","<<sv[9][0]<<","<<sv[47][0]<<endl;
            //sv[0][0] VIP 0,1,2 - saved in 2,3,4
            //sv[10][0] SOM 0,1,2 - saved in 5,6,7
            //sv[20][0] PV 0,1,2 - saved in 8
            //sv[30][0] ECS 0,1,...,9 - saved in 9,18
            //sv[40][0] PN 0,1,...,9 -saved in 19,28

            outputFile7 << times << "," << LFP[0] << "," << I[0] << "," << I[1] << "," << I[2] << endl;


            /********** PLASTICITY FOR ECS --> PN **********/
            int z;
            z = 0;
            if (stdp_rule == 4) {
                gAMPA_fromECStoPN[z] = gAMPA_fromECStoPN[z];
            }
            else
            {
                if (sv[40][z] > 0 && PN_svp[z] <= 0) // if post spikes
                {
                    gAMPA_fromECStoPN[z] = gAMPA_fromECStoPN[z] + sv[36][z];
                }
                if (sv[30][z] > 0 && ECS_svp[z] <= 0) // if pre spikes
                {
                    gAMPA_fromECStoPN[z] = gAMPA_fromECStoPN[z] + sv[46][z];
                }
            }
            if (gAMPA_fromECStoPN[z] < 0)
            {
                gAMPA_fromECStoPN[z] = 0;
            }
            if (gAMPA_fromECStoPN[z] > 0.185)
            {
                gAMPA_fromECStoPN[z] = 0.185;
            }

            /********** PLASTICITY FOR PN --> VIP **********/
            int zn; // PN index
            zn = 0;
            //z = 0;
            for (z = 0; z < VIPCOUNT; z++) { // z is the VIP index
                if (stdp_rule == 4) {
                    gAMPA_fromPNtoVIP[z][0] = gAMPA_fromPNtoVIP[z][0];
                } else {
                    if (sv[0][z] > 0 && VIP_svp[z] <= 0) // if post spikes
                    {
                        gAMPA_fromPNtoVIP[z][0] = gAMPA_fromPNtoVIP[z][0] + sv[47][z];
                    }
                    if (sv[40][0] > 0 && PN_svp[0] <= 0) // if pre spikes
                    {
                        gAMPA_fromPNtoVIP[z][0] = gAMPA_fromPNtoVIP[z][0] + sv[9][z];
                    }
                }
                if (gAMPA_fromPNtoVIP[z][0] < 0.01) {
                    gAMPA_fromPNtoVIP[z][0] = 0.01;
                }
                if (gAMPA_fromPNtoVIP[z][0] > 0.04) {
                    gAMPA_fromPNtoVIP[z][0] = 0.04;
                }
            }
            outputFile6<<gAMPA_fromECStoPN[0]<<','<<gAMPA_fromPNtoVIP[0][0]<<','<<gAMPA_fromPNtoVIP[1][0]<<','<<gAMPA_fromPNtoVIP[2][0]<<','<<gAMPA_fromPNtoVIP[0][1]<<endl;

            z = 0; // we consider only the 0th-PN and 0-th ECS
            if (sv[40][z] > 0 && PN_svp[z] <= 0)
            {
                sv[46][z] = sv[46][z] - A_minus; // M function
            }

            if (sv[30][z] > 0 && ECS_svp[z] <= 0) // if pre spikes
            {
                sv[36][z] = sv[36][z] + A_plus; // P function
            }

            for (z = 0; z < VIPCOUNT; z++) { // z is the VIP index
                if (sv[0][z] > 0 && VIP_svp[z] <= 0) // if VIP (post) fires
                {
                    sv[9][z] = sv[9][z] - A_minus_VIP; // M function
                }

                if (sv[40][0] > 0 && PN_svp[0] <= 0) // if PN (pre) spikes
                {
                    sv[47][z] = sv[47][z] + A_plus_VIP; // P function
                }
            }

            for (z=0; z<ECSCOUNT; z++){
                ECS_svp[z] = sv[30][z];
            }

            for (z=0; z<PNCOUNT; z++) {
                PN_svp[z] = sv[40][z];
            }

            for (z=0; z<PVCOUNT; z++) {
                PV_svp[z] = sv[20][z];
            }

            for (z=0; z<VIPCOUNT; z++) {
                VIP_svp[z] = sv[0][z];
            }



            itime2 = itime;

            times += TIMESTEP;

        }
        sleep(1);
    }
    return 0;
}


/*---------------------- Integration Methods ----------------------

 *-- IntRK4() -------------------------------------------------------
 *
 *--- general purpose Runge-Kutta integrator (fourth order)
 *    code is modified from "Numerical Recipes"
 *
 *--INPUTS
 *   int     svCount    = state variable count
 *   float  *stateVar   = state variable value array
 *   double *derivative = derivative value array
 *   float  *time      = time variable address
 *   float   dt         = timestep
 *-------------------------------------------------------------------
 */

int IntRK4( int svCount, int svCount2, double sv[SVCOUNT][SVCOUNT2], double times, double dt, double LFP[0], double I[CURRCOUNT], double gAMPA_fromECStoPN[ECSCOUNT], double gAMPA_fromPYRCStoPV, double isPoissonSpikeCS, double isPoissonSpikePN)

{
    int  i;                      /* loop counter */
    int  j;/* loop counter */

    /*-- temporary arrays --*/
    static double y [ SVCOUNT ][ SVCOUNT2 ];
    static double k1[ SVCOUNT ][ SVCOUNT2 ];
    static double k2[ SVCOUNT ][ SVCOUNT2 ];
    static double k3[ SVCOUNT ][ SVCOUNT2 ];
    static double k4[ SVCOUNT ][ SVCOUNT2 ];

    for ( i=0; i < SVCOUNT; i++ )       /* remember initial values */
    {
        for ( j=0; j < SVCOUNT2; j++)
        {
            y[i][j] = sv[i][j];
        }
    }

    /*-- pass 1 --*/
    Derive( SVCOUNT, SVCOUNT2, y, k1, times, LFP, I, gAMPA_fromECStoPN, gAMPA_fromPYRCStoPV, isPoissonSpikeCS,  isPoissonSpikePN);     /* get derivatives */

    for ( i=0; i < SVCOUNT; i++ )       /* update state variables */
    {
        for ( j=0; j < SVCOUNT2; j++)
        {
            y[ i ][j] = sv[ i ][j] + k1[ i ][j] * dt/2;
        }
    }

    /*-- pass 2 --*/
    times += dt/2;
    Derive( SVCOUNT, SVCOUNT2, y, k2, times, LFP, I, gAMPA_fromECStoPN, gAMPA_fromPYRCStoPV, isPoissonSpikeCS,  isPoissonSpikePN);     /* get derivatives */

    for ( i=0; i < SVCOUNT; i++ )       /* update state variables */
    {
        for ( j=0; j < SVCOUNT2; j++)
        {
            y[ i ][j] = sv[ i ][j] + k2[ i ][j] * dt/2;
        }
    }

    /*- pass 3 --*/
    Derive( SVCOUNT, SVCOUNT2, y, k3, times, LFP, I, gAMPA_fromECStoPN, gAMPA_fromPYRCStoPV, isPoissonSpikeCS,  isPoissonSpikePN);     /* get derivatives */

    for ( i=0; i < SVCOUNT; i++ )
    {
        for ( j=0; j < SVCOUNT2; j++)
        {
            y[ i ][j] = sv[ i ][j] + k3[ i ][j] * dt;    /* update state variables */
        }
    }

    /*-- pass 4 --*/
    times += dt/2;
    Derive( SVCOUNT, SVCOUNT2, y, k4, times, LFP, I, gAMPA_fromECStoPN, gAMPA_fromPYRCStoPV, isPoissonSpikeCS,  isPoissonSpikePN);  /* get derivatives */
    for ( i=0; i < SVCOUNT; i++ )
    {
        for ( j=0; j < SVCOUNT2; j++)
        {
            sv[ i ][j] = sv[ i ][j]  + ( k1[ i ][j] + 2*k2[ i ][j] + 2*k3[ i ][j] +k4[ i ][j] ) * dt/6;
        }
    }

    return 1;
}


void Derive( int svCount, int svCount2, double y[SVCOUNT][SVCOUNT2], double k[SVCOUNT][SVCOUNT2], double times, double LFP[0], double I[CURRCOUNT], double gAMPA_fromECStoPN[ECSCOUNT], double gAMPA_fromPYRCStoPV, double isPoissonSpikeCS, double isPoissonSpikePN) {
    int z;

    // synapses from VIP
    double GABA_gatingVariable_fromVIPtoSOM = 0;
    double GABA_gatingVariable_fromVIPtoPV = 0;
    for (z = 0; z < VIPCOUNT; z++) {
        GABA_gatingVariable_fromVIPtoSOM = GABA_gatingVariable_fromVIPtoSOM + y[7][z];
        GABA_gatingVariable_fromVIPtoPV = GABA_gatingVariable_fromVIPtoPV + y[8][z];
    }

    // synapses from SOM
    double GABA_gatingVariable_fromSOMtoPN = 0;
    double GABA_gatingVariable_fromSOMtoECS = 0;
    for (z = 0; z < SOMCOUNT; z++) {
        GABA_gatingVariable_fromSOMtoPN = GABA_gatingVariable_fromSOMtoPN + y[17][z];
        GABA_gatingVariable_fromSOMtoECS = GABA_gatingVariable_fromSOMtoECS + y[17][z];
    }

    double GABA_gatingVariable_fromPVtoPN = 0;
    double GABA_gatingVariable_fromPVtoECS = 0;
    double GABA_gatingVariable_fromPVtoSOM = 0;
    for (z = 0; z < PVCOUNT; z++) {
        GABA_gatingVariable_fromPVtoPN = GABA_gatingVariable_fromPVtoPN + y[24][z];
        GABA_gatingVariable_fromPVtoECS = GABA_gatingVariable_fromPVtoECS + y[24][z];
        GABA_gatingVariable_fromPVtoSOM = GABA_gatingVariable_fromPVtoSOM + y[24][z];
    }

    double GABA_gatingVariable_fromPNtoPV = 0;
    double AMPA_gatingVariable_fromPNtoVIP[VIPCOUNT];
    AMPA_gatingVariable_fromPNtoVIP[0] = 0;
    AMPA_gatingVariable_fromPNtoVIP[1] = 0;
    AMPA_gatingVariable_fromPNtoVIP[2] = 0;
    int zn;
    for (zn = 0; zn < PNCOUNT; zn++) {
        GABA_gatingVariable_fromPNtoPV = GABA_gatingVariable_fromPNtoPV + y[45][zn];
        for (z = 0; z < VIPCOUNT; z++) {
            AMPA_gatingVariable_fromPNtoVIP[z] = AMPA_gatingVariable_fromPNtoVIP[z] + gAMPA_fromPNtoVIP[z][zn] * y[45][zn];
        }
    }

    double AMPA_gatingVariable_fromECStoPN = 0;
    for (z = 0; z < ECSCOUNT; z++) {
        AMPA_gatingVariable_fromECStoPN = AMPA_gatingVariable_fromECStoPN + y[35][z];
    }

    double AMPA_gatingVariable_fromPYRCStoECS = 0;
    double AMPA_gatingVariable_fromPYRCStoPV = 0;
    for (z = 0; z < PYRCSCOUNT; z++) {
        AMPA_gatingVariable_fromPYRCStoECS = AMPA_gatingVariable_fromPYRCStoECS + y[55][z];
        AMPA_gatingVariable_fromPYRCStoPV = AMPA_gatingVariable_fromPYRCStoPV + y[55][z];
    }

    double AMPA_gatingVariable_fromPYRUStoPN = 0;
    double AMPA_gatingVariable_fromPYRUStoVIP = 0;
    for (z = 0; z < PYRUSCOUNT; z++) {
        AMPA_gatingVariable_fromPYRUStoPN = AMPA_gatingVariable_fromPYRUStoPN + y[65][z];
        AMPA_gatingVariable_fromPYRUStoVIP = AMPA_gatingVariable_fromPYRUStoVIP + y[65][z];
    }

    double gna = 112.5; // conductances in units of mS/cm^2
    double gk = 225;
    double gm = 0;

    double gd[VIPCOUNT];
    gd[0] = 3;
    gd[1] = 3;
    gd[2] = 3;

    double gl = 0.25;
    double Ena = 50;
    double Ek = -90;
    double El = -70;
    double Ei = -80;
    double Ee = 0;
    double Cm = 1;
    double taui = 5; // time constant of decay of GABAa current in units of ms

    /*** LFP that takes into account AMPA, D-, P-, and H- currents ****/
    LFP[0] = 0;
    I[0] = 0;
    I[1] = 0;
    I[2] = 0;

    /*** VIP equations ****/
    double ainfD[VIPCOUNT];
    double binfD[VIPCOUNT];
    double minf[VIPCOUNT];
    double hinf[VIPCOUNT];
    double tauh[VIPCOUNT];
    double ninf[VIPCOUNT];
    double taun[VIPCOUNT];
    double ahVIP[VIPCOUNT];
    double bhVIP[VIPCOUNT];

    double noise = 5 * sqrt(0.05) * box_muller(0, 1);
    for (z = 0; z < VIPCOUNT; z++) {
        if (condit != 1){
            noise = noise = 5 * sqrt(0.05) * box_muller(0, 1);
        }

        minf[z] = 1 / (1 + exp(-(y[0][z] + 24) / 11.5));

        k[0][z] = (1 / Cm) * (-gna * y[2][z] * (minf[z] * minf[z] * minf[z]) * (y[0][z] - Ena)
                              - gk * (y[3][z] * y[3][z]) * (y[0][z] - Ek)
                              - gd[z] * y[5][z] * y[5][z] * y[5][z] * y[6][z] * (y[0][z] - Ek)
                              - gl * (y[0][z] - El))
                  - AMPA_gatingVariable_fromPNtoVIP[z] * (y[0][z] - Ee)
                  + noise
                  + IappVIP[z];

        hinf[z] = 1 / (1 + exp((y[0][z] + 58.3) / 6.7));
        tauh[z] = 0.5 + 14 / (1 + exp((y[0][z] + 60) / 12));

        ninf[z] = 1 / (1 + exp(-(y[0][z] + 12.4) / 6.8));
        taun[z] = (0.087 + 11.4 / (1 + exp((y[0][z] + 14.6) / 8.6))) *
                  (0.087 + 11.4 / (1 + exp(-(y[0][z] - 1.3) / 18.7)));

        k[2][z] = (hinf[z] - y[2][z]) / tauh[z]; // derivative equation for h, Goloumb version
        k[3][z] = (ninf[z] - y[3][z]) / taun[z];

        ainfD[z] = 1 / (1 + exp(-(y[0][z] + 50) / 20));
        binfD[z] = 1 / (1 + exp((y[0][z] + 70) / 6));

        k[5][z] = (ainfD[z] - y[5][z]) / 2;
        k[6][z] = (binfD[z] - y[6][z]) / 150;


        /*** Synaptic equations ****/
        double sgi[VIPCOUNT];
        taui = 10;
        sgi[z] = 2 * (1 + tanh(y[0][z] / 4)); // i to i
        k[7][z] = sgi[z] * (1 - y[7][z]) - y[7][z] / taui;  // derivative equation for inhibitory GABAa synaptic activity (i to i kinetics)

        double sgiPV[VIPCOUNT];
        taui = 10;
        sgiPV[z] = 2 * (1 + tanh(y[0][z] / 4)); // i to i
        k[8][z] = sgiPV[z] * (1 - y[8][z]) - y[8][z] / taui;  // derivative equation for inhibitory GABAa synaptic activity (i to i kinetics)

        // M function - plasticity F -> VIP - VIP as post-synaptic neuron
        k[9][z] = - y[9][z]/tau_minus;
    }

    /*** SOM equations ****/

    double SOMgna = 52; // conductances in units of mS/cm^2
    double SOMgk = 11;
    double SOMgl = 0.62;
    double SOMgh[SOMCOUNT];
    SOMgh[0] = 1.5;
    SOMgh[1] = 1.4;
    SOMgh[2] = 1.45;

    double SOMgp = 0.5;
    double SOMEna = 55; // reversal potentials in units of mV
    double SOMEk = -90;
    double SOMEl = -65;
    double SOMEh = -20;

    double alpha_m[SOMCOUNT];
    double beta_m[SOMCOUNT];
    double alpha_n[SOMCOUNT];
    double beta_n[SOMCOUNT];
    double alpha_h[SOMCOUNT];
    double beta_h[SOMCOUNT];
    double h_f_inf[SOMCOUNT];
    double tauh_f[SOMCOUNT];
    double h_s_inf[SOMCOUNT];
    double tauh_s[SOMCOUNT];

    double IappSOM[SOMCOUNT];
    IappSOM[0] = 0.1;
    for (z = 1; z < SOMCOUNT; z++)
        IappSOM[z] = IappSOM[z - 1] + 0.00;

    for (z = 0; z < SOMCOUNT; z++) {
        k[10][z] = -SOMgna * y[12][z] * (y[11][z] * y[11][z] * y[11][z]) * (y[10][z] - SOMEna)
                   - SOMgk * (y[13][z] * y[13][z] * y[13][z] * y[13][z]) * (y[10][z] - SOMEk)
                   - SOMgl * (y[10][z] - SOMEl)
                   - SOMgh[z] * (0.65 * y[14][z] + 0.35 * y[15][z]) * (y[10][z] - SOMEh) // h current
                   - SOMgp * y[16][z] * (y[10][z] - SOMEna)// p current
                   - gGABA_fromVIPtoSOM * GABA_gatingVariable_fromVIPtoSOM * (y[10][z] - Ei)
                   + 4 * sqrt(0.05) * box_muller(0, 1)
                   + IappSOM[z];

        alpha_m[z] = -0.1 * (y[10][z] + 23) / (exp(-0.1 * (y[10][z] + 23)) - 1); //alpha_m * (1 - y) - beta_m * y
        beta_m[z] = 4 * exp(-(y[10][z] + 48) / 18);
        k[11][z] = alpha_m[z] * (1 - y[11][z]) - beta_m[z] * y[11][z];

        alpha_h[z] = 0.07 * exp(-(y[10][z] + 37) / 20); // alpha_h * (1 - y) - beta_h * y
        beta_h[z] = 1 / (exp(-0.1 * (y[10][z] + 7)) + 1);
        k[12][z] = alpha_h[z] * (1 - y[12][z]) - beta_h[z] * y[12][z];

        alpha_n[z] = -0.01 * (y[10][z] + 27) / (exp(-0.1 * (y[10][z] + 27)) - 1); // alpha_h * (1 - y) - beta_h * y
        beta_n[z] = 0.125 * exp(-(y[10][z] + 37) / 80);
        k[13][z] = alpha_n[z] * (1 - y[13][z]) - beta_n[z] * y[13][z];

        // H current
        h_f_inf[z] = 1 / (1 + exp((y[10][z] + 79.2) / 9.78));
        tauh_f[z] = 0.51 / (exp((y[10][z] - 1.7) / 10) + exp(-(y[10][z] + 340) / 52)) + 1;
        //(h_f_inf - y_O[4]) / tauh_f
        k[14][z] = (h_f_inf[z] - y[14][z]) / tauh_f[z];
        h_s_inf[z] = pow(1 / ((1 + exp((y[10][z] + 2.83) / 15.9))), 58);
        tauh_s[z] = 5.6 / (exp((y[10][z] - 1.7) / 14) + exp(-(y[10][z] + 260) / 43)) + 1;
        k[15][z] = (h_s_inf[z] - y[15][z]) / tauh_s[z];

        // P current
        k[16][z] = (1 / (1 + exp(-(y[10][z] + 38) / 6.5)) - y[16][z]) / 0.15; // (pinf - y_O[6]) / taup

        // Synaptic equations
        double sgiSOM[SOMCOUNT];
        sgiSOM[z] = 5 / 2 * (1 + tanh(y[10][z] / 0.1));
        k[17][z] = sgiSOM[z] * (1 - y[17][z]) - y[17][z] * 0.05;
    }

    /***PV equations ****/

    double PVgna = 100;
    double PVgk = 80;
    double PVgl = 0.1;
    double PVEna = 50;
    double PVEk = -100;
    double PVEl = -67;


    for (z = 0; z < PVCOUNT; z++) {
        double IappPV = 0;

        k[20][z] = -PVgna * y[22][z] * (y[21][z] * y[21][z] * y[21][z]) * (y[20][z] - PVEna)
                   - PVgk * (y[23][z] * y[23][z] * y[23][z] * y[23][z]) * (y[20][z] - PVEk)
                   - PVgl * (y[20][z] - PVEl)
                   - gGABA_fromVIPtoPV * GABA_gatingVariable_fromVIPtoPV * (y[20][z] - Ei)
                   - gAMPA_fromPNtoPV * GABA_gatingVariable_fromPNtoPV * (y[20][z] - Ee)
                   - gAMPA_fromPYRCStoPV * AMPA_gatingVariable_fromPYRCStoPV * (y[20][z] - Ee)
                   + 4 * sqrt(0.05) * box_muller(0, 1)
                   + IappPV;

        k[21][z] = 0.32 * (y[20][z] + 54) / (1 - exp(-(y[20][z] + 54) / 4)) * (1 - y[21][z]) -
                   0.28 * (y[20][z] + 27) / (exp((y[20][z] + 27) / 5) - 1) * y[21][z];     // derivative equation for m
        k[22][z] = 0.128 * exp(-(y[20][z] + 50) / 18) * (1 - y[22][z]) -
                   4 / (1 + exp(-(y[20][z] + 27) / 5)) * y[22][z];   // derivative equation for h
        k[23][z] = 0.032 * (y[20][z] + 52) / (1 - exp(-(y[20][z] + 52) / 5)) * (1 - y[23][z]) -
                   0.5 * exp(-(y[20][z] + 57) / 40) * y[23][z];  // derivative equation for n

        // Synaptic equations
        double sgiPV[PVCOUNT]; //double tauiPV = 16;
        sgiPV[z] = 15 / 2 * (1 + tanh(y[20][z] / 0.1));    //  2*(1+tanh(y[10][z]/4));
        k[24][z] = sgiPV[z] * (1 - y[24][z]) - y[24][z] * 0.12;
    }

    /*** ECS equations ****/

    double ECSgna = 100;
    double ECSgk = 80;
    double ECSgl = 0.1;
    double ECSEna = 50;
    double ECSEk = -100;
    double ECSEl = -67;
    double taue = 2;

    double alpha_m_ecs[ECSCOUNT];
    double beta_m_ecs[ECSCOUNT];
    double m_inf_ecs[ECSCOUNT];
    double alpha_h_ecs[ECSCOUNT];
    double beta_h_ecs[ECSCOUNT];
    double alpha_n_ecs[ECSCOUNT];
    double beta_n_ecs[ECSCOUNT];
    double tau_h_ecs[ECSCOUNT];
    double h_inf_ecs[ECSCOUNT];
    double tau_n_ecs[ECSCOUNT];
    double n_inf_ecs[ECSCOUNT];
    double w_inf_ecs[PNCOUNT];
    double tau_M_ecs[PNCOUNT];




    for (z = 0; z < ECSCOUNT; z++) {
        alpha_m_ecs[z] = 0.1 * (y[30][z] + 35) / (1 - exp(-(y[30][z] + 35) / 10));
        beta_m_ecs[z] = 4 * exp(-(y[30][z] + 60) / 18);
        m_inf_ecs[z] = alpha_m_ecs[z] / (alpha_m_ecs[z] + beta_m_ecs[z]);

        k[30][z] = -ECSgna * y[32][z] * (m_inf_ecs[z] * m_inf_ecs[z] * m_inf_ecs[z]) * (y[30][z] - ECSEna)
                   - ECSgk * (y[33][z] * y[33][z] * y[33][z] * y[33][z]) * (y[30][z] - ECSEk)
                   - ECSgl * (y[30][z] - ECSEl)
                   - gGABA_fromSOMtoECS * GABA_gatingVariable_fromSOMtoECS * (y[30][z] - Ei)
                   - gGABA_fromPVtoECS * GABA_gatingVariable_fromPVtoECS * (y[30][z] - Ei)
                   - gAMPA_fromPYRCStoECS[z] * AMPA_gatingVariable_fromPYRCStoECS * (y[30][z] - Ee)
                   + 4 * sqrt(0.05) * box_muller(0, 1)
                   + IappECS;

        //h
        alpha_h_ecs[z] = 0.07 * exp(-(y[30][z] + 58) / 20);
        beta_h_ecs[z] = 1. / (1 + exp(-(y[30][z] + 28) / 10));
        tau_h_ecs[z] = 1. / 5 * 1. / (alpha_h_ecs[z] + beta_h_ecs[z]);
        h_inf_ecs[z] = alpha_h_ecs[z] / (alpha_h_ecs[z] + beta_h_ecs[z]);
        k[32][z] = (h_inf_ecs[z] - y[32][z]) / tau_h_ecs[z];

        //n
        alpha_n_ecs[z] = -0.01 * (y[30][z] + 34) / (-1 + exp(-(y[30][z] + 34) / 10));
        beta_n_ecs[z] = 0.125 * exp(-(y[30][z] + 44) / 80);
        tau_n_ecs[z] = 1. / 5 * 1. / (alpha_n_ecs[z] + beta_n_ecs[z]);
        n_inf_ecs[z] = alpha_n_ecs[z] / (alpha_n_ecs[z] + beta_n_ecs[z]);
        k[33][z] = (n_inf_ecs[z] - y[33][z]) / tau_n_ecs[z];

        //M-current
        w_inf_ecs[z] = 1. / (1 + exp(-(y[30][z] + 35) / 10));
        tau_M_ecs[z] = 400. / (3.3 * exp((y[30][z] + 35) / 20) + exp(-(y[30][z] + 35) / 20));
        k[34][z] = (w_inf_ecs[z] - y[34][z]) / tau_M_ecs[z];

        // AMPA gating variable
        double sgeECS;
        sgeECS = 5 * (1 + tanh(y[30][z] / 4));
        k[35][z] = sgeECS * (1 - y[35][z]) - y[35][z] / taue;

        // P function - plasticity - pre-synaptic neuron
        k[36][z] = - y[36][z]/tau_plus;

    }

    /*** PN equations ****/
    double PNgna = 100;
    double PNgk = 80;
    double PNgl = 0.1;
    double PNEna = 50;
    double PNEk = -100;
    double PNEl = -67;

    double alpha_m_pn[PNCOUNT];
    double beta_m_pn[PNCOUNT];
    double m_inf_pn[PNCOUNT];
    double alpha_h_pn[PNCOUNT];
    double beta_h_pn[PNCOUNT];
    double alpha_n_pn[PNCOUNT];
    double beta_n_pn[PNCOUNT];
    double tau_h_pn[PNCOUNT];
    double h_inf_pn[PNCOUNT];
    double tau_n_pn[PNCOUNT];
    double n_inf_pn[PNCOUNT];
    double w_inf[PNCOUNT];
    double tau_M[PNCOUNT];


    for (z = 0; z < PNCOUNT; z++) {
        alpha_m_pn[z] = 0.1 * (y[40][z] + 35) / (1 - exp(-(y[40][z] + 35) / 10));
        beta_m_pn[z] = 4 * exp(-(y[40][z] + 60) / 18);
        m_inf_pn[z] = alpha_m_pn[z] / (alpha_m_pn[z] + beta_m_pn[z]);

        k[40][z] = -PNgna * y[42][z] * (m_inf_pn[z] * m_inf_pn[z] * m_inf_pn[z]) * (y[40][z] - PNEna)
                   - PNgk * (y[43][z] * y[43][z] * y[43][z] * y[43][z]) * (y[40][z] - PNEk)
                   - PNgl * (y[40][z] - PNEl)
                   - gGABA_fromPVtoPN * GABA_gatingVariable_fromPVtoPN * (y[40][z] - Ei)
                   - gGABA_fromSOMtoPN * GABA_gatingVariable_fromSOMtoPN * (y[40][z] - Ei)
                   - gAMPA_fromECStoPN[z] * AMPA_gatingVariable_fromECStoPN * (y[40][z] - Ee)
                   - gAMPA_fromPYRUStoPN[z] * AMPA_gatingVariable_fromPYRUStoPN *
                     (y[40][z] - Ee) // poisson spike train from PYRUS
                   + 4 * sqrt(0.05) * box_muller(0, 1)
                   + IappPN;

        //h
        alpha_h_pn[z] = 0.07 * exp(-(y[40][z] + 58) / 20);
        beta_h_pn[z] = 1. / (1 + exp(-(y[40][z] + 28) / 10));
        tau_h_pn[z] = 1. / 5 * 1. / (alpha_h_pn[z] + beta_h_pn[z]);
        h_inf_pn[z] = alpha_h_pn[z] / (alpha_h_pn[z] + beta_h_pn[z]);
        k[42][z] = (h_inf_pn[z] - y[42][z]) / tau_h_pn[z];

        //n
        alpha_n_pn[z] = -0.01 * (y[40][z] + 34) / (-1 + exp(-(y[40][z] + 34) / 10));
        beta_n_pn[z] = 0.125 * exp(-(y[40][z] + 44) / 80);
        tau_n_pn[z] = 1. / 5 * 1. / (alpha_n_pn[z] + beta_n_pn[z]);
        n_inf_pn[z] = alpha_n_pn[z] / (alpha_n_pn[z] + beta_n_pn[z]);
        k[43][z] = (n_inf_pn[z] - y[43][z]) / tau_n_pn[z];

        //M-current
        w_inf[z] = 1. / (1 + exp(-(y[40][z] + 35) / 10));
        tau_M[z] = 400. / (3.3 * exp((y[40][z] + 35) / 20) + exp(-(y[40][z] + 35) / 20));
        k[44][z] = (w_inf[z] - y[44][z]) / tau_M[z];

        // AMPA gating variable
        double sgePN;
        sgePN = 5 * (1 + tanh(y[40][z] / 4));
        k[45][z] = sgePN * (1 - y[45][z]) - y[45][z] / taue;

        // M function - plasticity - post-synaptic neuron
        k[46][z] = - y[46][z]/tau_minus;

        // P function - for plasticity F -> VIP - F pre-syn neuron
        k[47][z] = - y[47][z]/tau_plus;
    }

    for (z = 0; z < PYRCSCOUNT; z++) {
        alpha_m_pn[z] = 0.1 * (y[50][z] + 35) / (1 - exp(-(y[50][z] + 35) / 10));
        beta_m_pn[z] = 4 * exp(-(y[50][z] + 60) / 18);
        m_inf_pn[z] = alpha_m_pn[z] / (alpha_m_pn[z] + beta_m_pn[z]);

        k[50][z] = -PNgna * y[52][z] * (m_inf_pn[z] * m_inf_pn[z] * m_inf_pn[z]) * (y[50][z] - PNEna)
                   - PNgk * (y[53][z] * y[53][z] * y[53][z] * y[53][z]) * (y[50][z] - PNEk)
                   - PNgl * (y[50][z] - PNEl)
                   + 4 * sqrt(0.05) * box_muller(0, 1)
                   + 0.26
                   + strength_poisson_CS * isPoissonSpikeCS;

        //h
        alpha_h_pn[z] = 0.07 * exp(-(y[50][z] + 58) / 20);
        beta_h_pn[z] = 1. / (1 + exp(-(y[50][z] + 28) / 10));
        tau_h_pn[z] = 1. / 5 * 1. / (alpha_h_pn[z] + beta_h_pn[z]);
        h_inf_pn[z] = alpha_h_pn[z] / (alpha_h_pn[z] + beta_h_pn[z]);
        k[52][z] = (h_inf_pn[z] - y[52][z]) / tau_h_pn[z];

        //n
        alpha_n_pn[z] = -0.01 * (y[50][z] + 34) / (-1 + exp(-(y[50][z] + 34) / 10));
        beta_n_pn[z] = 0.125 * exp(-(y[50][z] + 44) / 80);
        tau_n_pn[z] = 1. / 5 * 1. / (alpha_n_pn[z] + beta_n_pn[z]);
        n_inf_pn[z] = alpha_n_pn[z] / (alpha_n_pn[z] + beta_n_pn[z]);
        k[53][z] = (n_inf_pn[z] - y[53][z]) / tau_n_pn[z];

        //M-current
        w_inf[z] = 1. / (1 + exp(-(y[50][z] + 35) / 10));
        tau_M[z] = 400. / (3.3 * exp((y[50][z] + 35) / 20) + exp(-(y[50][z] + 35) / 20));
        k[54][z] = (w_inf[z] - y[54][z]) / tau_M[z];

        // AMPA gating variable
        double sgeECS;
        sgeECS = 5 * (1 + tanh(y[50][z] / 4));
        k[55][z] = sgeECS * (1 - y[55][z]) - y[55][z] / taue;

    }

    /***PYR US equations ****/
    for (z = 0; z < PYRUSCOUNT; z++) {
        alpha_m_pn[z] = 0.1 * (y[60][z] + 35) / (1 - exp(-(y[60][z] + 35) / 10));
        beta_m_pn[z] = 4 * exp(-(y[60][z] + 60) / 18);
        m_inf_pn[z] = alpha_m_pn[z] / (alpha_m_pn[z] + beta_m_pn[z]);

        k[60][z] = -PNgna * y[62][z] * (m_inf_pn[z] * m_inf_pn[z] * m_inf_pn[z]) * (y[60][z] - PNEna)
                   - PNgk * (y[63][z] * y[63][z] * y[63][z] * y[63][z]) * (y[60][z] - PNEk)
                   - PNgl * (y[60][z] - PNEl)
                   + 4 * sqrt(0.05) * box_muller(0, 1)
                   + 0.26
                   + strength_poisson_PN * isPoissonSpikePN;

        //h
        alpha_h_pn[z] = 0.07 * exp(-(y[60][z] + 58) / 20);
        beta_h_pn[z] = 1. / (1 + exp(-(y[60][z] + 28) / 10));
        tau_h_pn[z] = 1. / 5 * 1. / (alpha_h_pn[z] + beta_h_pn[z]);
        h_inf_pn[z] = alpha_h_pn[z] / (alpha_h_pn[z] + beta_h_pn[z]);
        k[62][z] = (h_inf_pn[z] - y[62][z]) / tau_h_pn[z];

        //n
        alpha_n_pn[z] = -0.01 * (y[60][z] + 34) / (-1 + exp(-(y[60][z] + 34) / 10));
        beta_n_pn[z] = 0.125 * exp(-(y[60][z] + 44) / 80);
        tau_n_pn[z] = 1. / 5 * 1. / (alpha_n_pn[z] + beta_n_pn[z]);
        n_inf_pn[z] = alpha_n_pn[z] / (alpha_n_pn[z] + beta_n_pn[z]);
        k[63][z] = (n_inf_pn[z] - y[63][z]) / tau_n_pn[z];

        //M-current
        w_inf[z] = 1. / (1 + exp(-(y[60][z] + 35) / 10));
        tau_M[z] = 400. / (3.3 * exp((y[60][z] + 35) / 20) + exp(-(y[60][z] + 35) / 20));
        k[64][z] = (w_inf[z] - y[64][z]) / tau_M[z];


        // AMPA gating variable
        double sgePN;
        sgePN = 5 * (1 + tanh(y[60][z] / 4));
        k[65][z] = sgePN * (1 - y[65][z]) - y[65][z] / taue;
    }


    for (z = 0; z < VIPCOUNT; z++) // counting each VIP
    {
        LFP[0] = LFP[0] - AMPA_gatingVariable_fromPNtoVIP[z] * (y[0][z] - Ee);
    }

    for (z = 0; z < PNCOUNT; z++) // counting each PN
    {
        LFP[0] = LFP[0] - gAMPA_fromECStoPN[z] * AMPA_gatingVariable_fromECStoPN * (y[40][z] - Ee);
    }

    for (z = 0; z < PVCOUNT; z++)  // counting each PV
    {
        LFP[0] = LFP[0] - gAMPA_fromPNtoPV * GABA_gatingVariable_fromPNtoPV * (y[20][z]-Ee);
    }

    for (z = 0; z < SOMCOUNT; z++){
        I[0] = I[0] - SOMgh[z] * (0.65 * y[14][z] + 0.35 * y[15][z]) * (y[10][z] - SOMEh); // h current
        I[1] = I[1] - SOMgp * y[16][z] * (y[10][z] - SOMEna);// p current
    }

    for (z = 0; z < VIPCOUNT; z++) {
        I[2] = I[2] -gd[z] * y[5][z] * y[5][z] * y[5][z] * y[6][z] * (y[0][z] - Ek); // d current
    }




}


double rnd()  {   // returns a uniform random number between 0 and 1
    return((double)rand()/RAND_MAX);
}


double randgauss( double min, double max, double sigma, double centre)
{
    double gauss = (min + (max-min) * (double)rand()/RAND_MAX); //create random domain between [min,max]
    return gauss;
}



float box_muller(float m, float s)    /* normal random variate generator */  //generates Gaussian random numbers
{                        /* mean m, standard deviation s */
    float x1, x2, w, y1;
    static float y2;
    static int use_last = 0;

    if (use_last)                /* use value from previous call */
    {
        y1 = y2;
        use_last = 0;
    }
    else
    {
        do {
            x1 = 2.0 * randgauss(0,1,1,0) - 1.0;
            x2 = 2.0 * randgauss(0,1,1,0) - 1.0;
            w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 );

        w = sqrt( (-2.0 * log( w ) ) / w );
        y1 = x1 * w;
        y2 = x2 * w;
        use_last = 1;
    }

    return( m + y1 * s );

}
