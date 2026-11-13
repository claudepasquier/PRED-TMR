/*
********************************************************************************
* PRED-TMR.cpp                                                 Claude Pasquier *
* Created 10/06/98                                        University of Athens *
*                                                                              *
* Predict transmembrane domains in a protein sequence                          *
*                                                                              *
*                                                                              *
*                                                                              *
********************************************************************************
*/

#include "cgipp.h"
#include "PRED-TMR.h"
#include "writeLog.h"

#include <string.h>
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>

const char AA[22] = "ACDEFGHIKLMNPQRSTVWYX";

// values computed from Swissprot 35 (69113 sequences)
const double freqOfAA[20] =     {0.0757, 0.0168, 0.0529, 0.0634, 0.0408, 0.0683, 0.0223, 0.0578, 0.0595, 0.0939, 0.0235, 0.0450, 0.0491, 0.0400, 0.0515, 0.0715, 0.0570, 0.0655, 0.0124, 0.0318};

// value when we proceed all the TM sequences in SP35 (excluding non ended TM)
//const double freqOfAA[20] =     {0.0752, 0.0187, 0.0433, 0.0505, 0.0499, 0.0689, 0.0204, 0.0654, 0.0468, 0.1058, 0.0257, 0.0430, 0.0478, 0.0353, 0.0452, 0.0742, 0.0605, 0.0711, 0.0171, 0.0341};

// value when we proceed all the TM sequences in SP35 (excluding non ended TM)
const double freqOfAAInTM[20] = {0.1047, 0.0202, 0.0081, 0.0083, 0.0912, 0.0791, 0.0088, 0.1204, 0.0064, 0.1732, 0.0353, 0.0175, 0.0293, 0.0109, 0.0064, 0.0576, 0.0501, 0.1150, 0.0222, 0.0342};

float table1[20][10] = {
    1.06902,   0.989605,    1.04213,   0.922376,   0.793477,    1.30621,    1.43313,    1.36232,    1.49283,    1.52134,   
   0.617683,   0.657246,   0.698084,   0.787419,   0.740837,    0.63555,   0.791247,   0.807838,   0.873563,    1.05287,   
   0.749396,   0.664802,   0.628002,   0.795278,     1.1542,   0.101321,   0.138839,   0.215786,   0.168232,   0.154133,   
    0.69983,   0.639244,   0.610494,   0.644839,    0.92809,   0.0984046,  0.146449,   0.185811,   0.125418,   0.140082,   
    1.21061,    1.35475,    1.24021,   0.993785,   0.396258,    2.59287,    2.33389,    2.16701,    2.20977,    2.05875,   
     1.0399,    1.03483,   0.976442,    1.04985,   0.865299,   0.937955,    1.10824,     1.1364,    1.22088,    1.13827,   
   0.992616,   0.996372,   0.993689,    1.09724,    1.42722,   0.300467,   0.379877,   0.628835,   0.556401,   0.442653,   
   0.996274,    1.05736,     1.0528,   0.827607,   0.301093,    2.71576,    2.18674,    2.11153,    2.10811,    2.19085,   
   0.999763,   0.970133,   0.941309,    1.11062,    1.59539,   0.086673,   0.123156,   0.173749,   0.116908,   0.118923,   
    1.01348,    1.06762,    1.09744,   0.891112,   0.534479,    2.02294,    1.90231,     1.7696,    1.66832,     1.9093,   
    1.42454,    1.08239,      1.091,    1.12543,   0.729484,     1.6408,    1.58593,    1.81188,    1.53805,    1.48856,   
   0.937446,   0.986294,   0.996174,    1.11994,    1.44267,   0.276898,   0.374594,   0.426187,   0.352914,   0.417679,   
   0.925405,   0.902998,   0.914081,   0.909022,   0.936969,   0.514622,   0.514863,   0.601597,   0.600634,   0.563049,   
   0.901526,   0.902744,   0.882041,   0.890566,      1.213,   0.178418,   0.257275,   0.362316,   0.268844,   0.319081,   
    1.31799,    1.28724,    1.29549,    1.69623,    2.38123,   0.153249,   0.125999,   0.170248,   0.157748,   0.163248,   
   0.987555,   0.992227,   0.967736,    1.00898,    1.19654,   0.409588,   0.700746,   0.716698,   0.780183,   0.720726,   
    1.00485,    1.10271,     1.0757,    1.14897,    1.22712,    0.65175,    0.80384,   0.862949,   0.944639,     1.0053,   
   0.931682,    0.93327,   0.959481,   0.788117,   0.419177,    1.99204,    1.75872,    1.66202,    1.75157,    1.78374,   
    1.66315,    1.84874,    1.77081,    1.83336,    1.04485,    2.87411,    2.51318,    2.10098,    2.09176,    1.45397,   
    1.14929,    1.17764,    1.46939,     1.3367,    1.14338,    1.37371,    1.42096,    1.37568,    1.45718,    1.02841   
};
float table2[20][10] = {
    1.52134,    1.49283,    1.36232,    1.43313,    1.30621,   0.793477,   0.922376,    1.04213,   0.989605,    1.06902,   
    1.05287,   0.873563,   0.807838,   0.791247,    0.63555,   0.740837,   0.787419,   0.698084,   0.657246,   0.617683,   
   0.154133,   0.168232,   0.215786,   0.138839,   0.101321,     1.1542,   0.795278,   0.628002,   0.664802,   0.749396,   
   0.140082,   0.125418,   0.185811,   0.146449,   0.0984046,   0.92809,   0.644839,   0.610494,   0.639244,    0.69983,   
    2.05875,    2.20977,    2.16701,    2.33389,    2.59287,   0.396258,   0.993785,    1.24021,    1.35475,    1.21061,   
    1.13827,    1.22088,     1.1364,    1.10824,   0.937955,   0.865299,    1.04985,   0.976442,    1.03483,     1.0399,   
   0.442653,   0.556401,   0.628835,   0.379877,   0.300467,    1.42722,    1.09724,   0.993689,   0.996372,   0.992616,   
    2.19085,    2.10811,    2.11153,    2.18674,    2.71576,   0.301093,   0.827607,     1.0528,    1.05736,   0.996274,   
   0.118923,   0.116908,   0.173749,   0.123156,   0.086673,    1.59539,    1.11062,   0.941309,   0.970133,   0.999763,   
     1.9093,    1.66832,     1.7696,    1.90231,    2.02294,   0.534479,   0.891112,    1.09744,    1.06762,    1.01348,   
    1.48856,    1.53805,    1.81188,    1.58593,     1.6408,   0.729484,    1.12543,      1.091,    1.08239,    1.42454,   
   0.417679,   0.352914,   0.426187,   0.374594,   0.276898,    1.44267,    1.11994,   0.996174,   0.986294,   0.937446,   
   0.563049,   0.600634,   0.601597,   0.514863,   0.514622,   0.936969,   0.909022,   0.914081,   0.902998,   0.925405,   
   0.319081,   0.268844,   0.362316,   0.257275,   0.178418,      1.213,   0.890566,   0.882041,   0.902744,   0.901526,   
   0.163248,   0.157748,   0.170248,   0.125999,   0.153249,    2.38123,    1.69623,    1.29549,    1.28724,    1.31799,   
   0.720726,   0.780183,   0.716698,   0.700746,   0.409588,    1.19654,    1.00898,   0.967736,   0.992227,   0.987555,   
     1.0053,   0.944639,   0.862949,    0.80384,    0.65175,    1.22712,    1.14897,     1.0757,    1.10271,    1.00485,   
    1.78374,    1.75157,    1.66202,    1.75872,    1.99204,   0.419177,   0.788117,   0.959481,    0.93327,   0.931682,   
    1.45397,    2.09176,    2.10098,    2.51318,    2.87411,    1.04485,    1.83336,    1.77081,    1.84874,    1.66315,   
    1.02841,    1.45718,    1.37568,    1.42096,    1.37371,    1.14338,     1.3367,    1.46939,    1.17764,    1.14929   

};

int nbInputFile = 0;
int fileName[100] = {0};
char outputFormat = 'C';
char outputMode ='H'; // 'H' == HTML mode, 'T' = Text mode
int hydroCutoff = 7.5;
bool showObservedTM = FALSE;
bool includeSPInfo = TRUE;
bool call_CoPreTHi = FALSE;
bool clearForm = FALSE;
char seqName[20];
char* seqData;

int TMInSwiss[10000] = {0};
int TM[10000] = {0};
float startInd[10000] = {0};
float endInd[10000] = {0};
float hydroInd[10000] = {0};
int nbAA;

bool main(int argc, char *argv[])
{
    bool error = FALSE;
    if (getenv("HOME")) { // the program is called locally
        if (!ParseArg(argc, argv)) return FALSE;
        outputMode = 'T';
        for (int i = 0 ; i < nbInputFile ; i++ ) {
            error |= ProcessFile(argv[fileName[i]]);
        }
    }
    else { // the program is called by a web browser with a method 'GET' or 'POST'
        writeLog();
        if (!ParseArg()) return FALSE;
        cout << "Content-type: text/html" << endl << endl;
        ProcessFile("tempfile");
        system("rm -f tempfile");
    }
    return error;
}


bool ProcessFile(char* fileName) {
    /* open the file */

    SeqFile fileIn = SeqFile();
    fileIn.open(fileName, SeqFile::seek);
    if (!fileIn.good()) {
        cerr << "PRED-TMR: cannot open the file '" << fileName << "'" << endl;
        return FALSE;
    }
    if ((outputMode == 'H') && (fileIn.countSequences() > 1)) {
        cout << "PRED-TMR: ***ERROR***<BR>Only one sequence can be submitted at a time" << endl;
        return FALSE;
    }

    if (fileIn.getFormat() != SeqFile::swiss) showObservedTM = FALSE;
    Sequence* seq;
    int count = 0;
    while ((seq = fileIn.getSequence(count))) {
        for (int i = 0 ; i < 10000 ; i++ ) {
            TMInSwiss[i] = 0;
            TM[i]        = 0;
            startInd[i]  = 0;
            endInd[i]    = 0;
            hydroInd[i]  = 0;
        }


        nbAA = strlen(seq->GetSequenceData());
        seqData = new char[nbAA+1];
        strcpy(seqData, seq->GetSequenceData());
        if (strlen(seq->GetEntryName()) > 0) strcpy(seqName, seq->GetEntryName());

        if (strlen(seqData) < 11) {
            DisplayResult();
            count++;
            break;
        }
        for (i = 0 ; i < seq->GetNbTransmem() ; i++ ) {
            int from;
            int to;
            int extent;
            int dir;
            seq->GetTransmem(i+1, from, to, extent, dir);
            for (int j = from-1 ; j <= to-1 ; j++ ) TMInSwiss[j] = 1;
        }


        ComputeIndicators();

        MarkTM();

        if (call_CoPreTHi) Generate_CoPreTHi_Page();
        else DisplayResult();

        delete[] seqData;
        delete seq;
        count++;
    }
    fileIn.close();
    return TRUE;
}

bool ComputeIndicators() {
        for (int i = 0 ; i < 5 ; i++ ) {
            hydroInd[i] = freqOfAAInTM[strchr(AA, seqData[i])-AA]/freqOfAA[strchr(AA, seqData[i])-AA];
        }
        for (i = 0 ; i < strlen(seqData) - 10 ; i++ ) {
            float indice1 = 1;
            for (int j = 0 ; j < 10 ; j++ ) {
                indice1 *= table1[strchr(AA, seqData[i+j])-AA][j];
            }
            float indice2 = 1;
            for (j = 0 ; j < 10 ; j++ ) {
                indice2 *= table2[strchr(AA, seqData[i+j])-AA][j];
            }
            float indice3 = 1;
            indice3 = freqOfAAInTM[strchr(AA, seqData[i+5])-AA]/freqOfAA[strchr(AA, seqData[i+5])-AA];

            startInd[i+5] = indice1;
            endInd[i+5] = indice2;
            hydroInd[i+5] = indice3;
        }
        for (i = strlen(seqData)-5 ; i < strlen(seqData) ; i++ ) {
            hydroInd[i] = freqOfAAInTM[strchr(AA, seqData[i])-AA]/freqOfAA[strchr(AA, seqData[i])-AA];
       }
    return TRUE;
}

bool MarkTM() {
    int bestStr[10000] = {0};
    int endIndx[10000] = {0};
//    int bestStart[10000] = {0};
    float bestMiddle[10000] = {0};
//    int bestEnd[10000] = {0};
//    int bestLen[10000] = {0};

    for (int i = 0 ; i < nbAA ; i++ ) {
        float startStr = 0;
        float middleStr = 0;
        float totalTM = 0;
        int nbAAInTM = 0;
        // test if it is a good start
        startStr = fsqrt(startInd[i]/endInd[i]);
        if (startStr>startInd[i]) startStr = startInd[i];
        startStr *= startInd[i];
        startStr = sqrtf(startStr);

        for (int j = i ; j < i+40 ; j++ ) {
            middleStr += log(hydroInd[j]);
            totalTM += hydroInd[j];
            nbAAInTM++;
            float endStr = fsqrt(endInd[j+1]/startInd[j+1]);
            if (endStr>endInd[j+1]) endStr = endInd[j+1];
            endStr *= endInd[j+1];
            endStr = fsqrt(endStr);
            int length = j-i+1;
            float normMiddle = fexp(middleStr*10/length);
            float lengthPenalty = fexp(abs(j-i+1-21));
            float totStr = startStr + normMiddle + endStr - lengthPenalty;
            if (startStr < 1) continue;
            if (endStr < 1) continue;
            if (normMiddle < 1) continue;
            if (totStr >  bestStr[i]) {
                bestStr[i] = totStr;
                endIndx[i] = j;
//                bestStart[i] = startStr;
                bestMiddle[i] = normMiddle;
//                bestEnd[i] = endStr;
//                bestLen[i] = lengthPenalty;
            }
        }
    }

    // select the TM part
    
    int imax = -1;
    while (1) {
        imax = -1;
        for (int i = 0 ; i < nbAA ; i++ ) {
            if (TM[i] > 0) continue;
//            if (bestStr[i] <10) continue;
            if (bestMiddle[i]<hydroCutoff) continue;
//            if ((bestMiddle[i]<7.5) && (bestStr[i] < 30)) continue;
//            if (bestMiddle[i]<3.5) continue;
            // get the TM mark
            if ((imax == -1) || (bestStr[i]>bestStr[imax])) {
                if (TM[endIndx[i]] == 0) imax = i;
            }
        }
        if (imax == -1) break;
        for (i = imax ; i <= endIndx[imax] ; i++ ) {
                TM[i] = 1;
            }
    }
    return TRUE;
}

// Display result
void DisplayResult() {
    if (outputFormat == 'T') {
        if (strlen(seqName) > 0) cout << setw(10) << seqName << " | ";

        // count the number of predicted transmembrane regions
        int nbTM = 0;
        if (TM[0] == 1 ) nbTM++;
        for (int i = 1 ; i < nbAA ; i++ ) {
            if ( (TM[i-1] == 0) && (TM[i] == 1) ) nbTM++;
        }
        cout << setw(2) << nbTM << " | ";

        // count the number of observed transmembrane regions
        int nbObsTM = 0;
        if (TMInSwiss[0] == 1 ) nbObsTM++;
        for (i = 1 ; i < nbAA ; i++ ) {
            if ( (TMInSwiss[i-1] == 0) && (TMInSwiss[i] == 1) ) nbObsTM++;
        }
        cout << setw(2) << nbObsTM << " | ";

        // print the rating
        float Ca;
        float Qa;
        float SM;
        if (rating(Qa, Ca, SM)) {
            cout << setprecision(3) << Qa << " | ";
            cout << setprecision(3) << Ca << " | ";
            cout << setw(2) << SM << " | ";
        }
        else cout << "   |   |   | ";

        // print the detail of predicted transmembrane regions
        int inTM = 0;
        for (i = 0 ; i < nbAA ; i++ ) {
            if ((inTM==0) && (TM[i]==1)) {
                cout << "(" << i+1 << "-";
                inTM = 1;
                nbTM--;
            }
            else if ((inTM==1) && (TM[i+1]==0)) {
                cout << i+1 << ") ";
                inTM = 0;
            }
        }
        cout << " | ";

        // print the detail of observed transmembrane regions
        inTM = 0;
        for (i = 0 ; i < nbAA ; i++ ) {
            if ((inTM==0) && (TMInSwiss[i]==1)) {
                cout << "(" << i+1 << "-";
                inTM = 1;
                nbObsTM--;
            }
            else if ((inTM==1) && (TMInSwiss[i+1]==0)) {
                cout << i+1 << ") ";
                inTM = 0;
            }
        }
        cout << " | " << endl;
        return;
    }
    if (outputMode == 'H') {
        cout << "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
        cout << "<HTML>" << endl;
        cout << "  <HEAD BASE href=\"http://o2.db.uoa.gr/PRED-TMR/index.html\">" << endl;
        cout << "    <TITLE>PRED-TMR: Prediction of Transmembrane regions in proteins</TITLE>" << endl;
        cout << "  </HEAD>" << endl;
        cout << "  <BODY BGCOLOR=\"#FFF0C0\">" << endl;
        cout << "    <H2>PRED-TMR: Prediction of Transmembrane regions in proteins</H2>" << endl;
        cout << "    <P>" << endl;
        cout << "    <P>" << endl;
        cout << "    <H3>Sequence: ";
        if (strlen(seqName) > 0) cout << seqName;
        cout << "</H3>" << endl;
        cout << "    <PRE>";
    }
    else {
        cout << "PRED-TMR: Prediction of Transmembrane regions in proteins" << endl;
        cout << endl;
        cout << "Sequence: ";
        if (strlen(seqName) > 0) cout << seqName;
        cout << endl << endl;
    }

    int inTM = 0;
    int inObsTM;
    cout << "               1         2         3         4         5         6" << endl;
    cout << "      123456789012345678901234567890123456789012345678901234567890" << endl;
    cout << "          .    |    .    |    .    |    .    |    .    |    .    |" << endl;
    for ( int line = 0 ; line*60 <= strlen(seqData) ; line++) {
        if (showObservedTM) {
           cout << "      ";
            if ((outputMode == 'H') && (inObsTM==1)) cout << "<FONT COLOR=RED><B>";
            for (int res = 0 ; (line*60+res < strlen(seqData)) && (res<60) ; res++) {
                int indx = line*60+res;
                if ((inObsTM==0) && (TMInSwiss[indx]==1)) {
                    if (outputMode == 'H') cout << "<FONT COLOR=RED><B>";
                    inObsTM = 1;
                }
                else if ((inObsTM==1) && (TMInSwiss[indx]==0)) {
                    if (outputMode == 'H') cout << "</B></FONT>";
                    inObsTM = 0;
                }
                if (TMInSwiss[indx]==1)
                    cout << '_';
                else cout << ' ';
            }
            if (inObsTM==1) cout << "</B></FONT>";
            cout << endl;
        }
        cout << " " << setw(4) << setfill('0') << line*60 << " " << setfill(' ');
        if ((outputMode == 'H') && (inTM==1)) cout << "<FONT COLOR=BLUE><B>";
        for (int res = 0 ; (line*60+res < strlen(seqData)) && (res<60) ; res++) {
            int indx = line*60+res;
            if ((inTM==0) && (TM[indx]==1)) {
                if (outputMode == 'H') cout << "<FONT COLOR=BLUE><B>";
                inTM = 1;
            }
            else if ((inTM==1) && (TM[indx]==0)) {
                if (outputMode == 'H') cout << "</B></FONT>";
                inTM = 0;
            }
            cout << seqData[indx];
        }
        if ((outputMode == 'H') && (inTM==1)) cout << "</B></FONT>";
        cout << endl;
    }
    if (outputMode == 'H') {
        cout << endl;
        cout << "<I>(predicted transmembrane regions are indicated in blue)</I>" << endl;
        if (showObservedTM) cout << "<I>(observed transmembrane regions are overlined in red)</I>" << endl;
    }
    cout << endl << endl;

    if (showObservedTM) {
        float Ca;
        float Qa;
        float SM;
        if (rating(Qa, Ca, SM)) {
            if (outputMode == 'H') {
                cout << "<H3>Rating of the prediction against the observed structure</H3>" << endl;
                cout << "Agreement factor (Q)        = <FONT SIZE=+2><B>" << setprecision(4) << Qa << "%</B></FONT>" << endl << endl;
                char sep[3] = "  ";
                if (Qa < 100) sep[1] = '\x0';
                if (Qa < 10) sep[0] = '\x0';
                cout << "Correlation coefficient (C) = <FONT SIZE=+2><B>" << sep << setprecision(3) << Ca << "</B></FONT>" << endl << endl;
                cout << "Ratio of segment matches     = <FONT SIZE=+2><B>" << sep << setprecision(3) << SM << "</B></FONT>" << endl << endl;
            }
            else {
                cout << "Rating of the prediction against the observed structure" << endl << endl;
                cout << "Agreement factor (Q)        = " << setprecision(4) << Qa << "%" << endl;
                char sep[3] = "  ";
                if (Qa < 100) sep[1] = '\x0';
                if (Qa < 10) sep[0] = '\x0';
                cout << "Correlation coefficient (C) = " << sep << setprecision(3) << Ca << endl;
                cout << "Ratio of segment matches     = " << sep << setprecision(3) << SM << endl;
            }
        }
        else {
            if (outputMode == 'H') {
                cout << "<H3>Rating of the prediction against the observed structure</H3>" << endl;
                cout << "Agreement factor (Q)        = <FONT SIZE=+2><B>na</B></FONT>" << endl << endl;
                cout << "Correlation coefficient (C) = <FONT SIZE=+2><B>na</B></FONT>" << endl << endl;
                cout << "Ratio of segment matches     = <FONT SIZE=+2><B>na</B></FONT>" << endl << endl;
            }
            else {
                cout << "Rating of the prediction against the observed structure" << endl << endl;
                cout << "Agreement factor (Q)        = na" << endl;
                cout << "Correlation coefficient (C) = na" << endl;
                cout << "Ratio of segment matches     = na" << endl;
            }
        }
    }

    if (outputMode == 'H') {
        cout << "<H3>Detail of predicted transmembrane segments:</H3>" << endl;
    }
    else {
        cout << endl <<endl << "Detail of predicted transmembrane segments:" << endl;
    }


    if (outputFormat == 'C') {
        int nbObsTM = 0;
        for (int i = 1 ; i < nbAA ; i++ ) {
            if ( (TMInSwiss[i-1] == 0) && (TMInSwiss[i] == 1) ) nbObsTM++;
        }
        int nbTM = 0;
        for (i = 1 ; i < nbAA ; i++ ) {
            if ( (TM[i-1] == 0) && (TM[i] == 1) ) nbTM++;
        }
        cout << "Number of predicted transmembrane segments: ";
        if  (outputMode == 'H') cout << "<FONT SIZE=+2><B>";
        cout << setw(2) << nbTM;
        if  (outputMode == 'H') cout << "</B></FONT>";
        cout << endl << endl;
        if (showObservedTM) {
            cout << "Number of observed transmembrane segments:  ";
            if  (outputMode == 'H') cout << "<FONT SIZE=+2><B>";
            cout << setw(2) << nbObsTM;
            if  (outputMode == 'H') cout << "</B></FONT>";
            cout << endl << endl;
            if ((nbObsTM > 0) || (nbTM > 0)) {
                if  (outputMode == 'H') cout << "<B>";
                cout << "     predicted            observed " << endl;
                cout << "#TM  begin  end          begin  end" << endl;
                if  (outputMode == 'H') cout << "</B>";
                nbTM = 0;
                nbObsTM = 0;
                for (int i = 0 ; i < nbAA ; i++ ) {
                    if (TM[i]==1) { // beginning of predicted TM region
                        for (int j = i; j < nbAA ; j++ ) {
                            if (TM[j]==0) break;
                            TM[j] = 0;
                        }
                        cout << setw(2) << ++nbTM << "   " << setw(4) << i+1 << "  " << setw(4) << j;
                        bool match = FALSE;
                        for (int k=i ; k < j ; k++ ) {
                            if (TMInSwiss[k] == 1) {
                                match = TRUE;
                                break;
                            }
                        }
                        for (int l = k; l < nbAA ; l++ ) {
                            if (TMInSwiss[l]==0) break;
                            TMInSwiss[l] = 0;
                        }
                        if (match) cout << "          " << setw(4) << k+1 << "  " << setw(4) << l;
                        cout << endl;
                    }
                    else if (TMInSwiss[i]==1) { // beginning of observed TM region
                        for (int j = i; j < nbAA ; j++ ) {
                            if (TMInSwiss[j]==0) break;
                            TMInSwiss[j] = 0;
                        }
                        int match = FALSE;
                        for (int k=i ; k < j ; k++ ) {
                            if (TM[k] == 1) {
                                match = TRUE;
                                break;
                            }
                        }
                        for (int l = k; l < nbAA ; l++ ) {
                            if (TM[l]==0) break;
                            TM[l] = 0;
                        }
                        if (match) {
                            cout << setw(2) << ++nbTM << "   " << setw(4) << k+1 << "  " << setw(4) << l;
                        }
                        else cout << "               ";
                        cout << "          " << setw(4) << i+1 << "  " << setw(4) << j << endl;
                    }
                }
            }
        }
        else if (nbTM > 0) {
            nbTM = 0;
            int inTM = 0;
            if  (outputMode == 'H') cout << "<B>";
            cout << "#TM  begin  end" << endl;
            if  (outputMode == 'H') cout << "</B>";
            for (int i = 0 ; i < nbAA ; i++ ) {
                if ((inTM==0) && (TM[i]==1)) {
                    cout << setw(2) << ++nbTM << "   " << setw(4) << i+1 << "  ";
                    inTM = 1;
                }
                else if ((inTM==1) && (TM[i+1]==0)) {
                    cout << setw(4) << i+1 << endl;
                    inTM = 0;
                }
            }
        }
    }
    else if (outputFormat == 'R') {
        cout << " ";
        int inTM = 0;
        int nbTM = 0;
        // count the number of transmembrane regions
        if (TM[0] == 1 ) nbTM++;
        for (int i = 1 ; i < nbAA ; i++ ) {
            if ( (TM[i-1] == 0) && (TM[i] == 1) ) nbTM++;
        }
        cout << setw(2) << nbTM << " ";
        for (i = 0 ; i < nbAA ; i++ ) {
            if ((inTM==0) && (TM[i]==1)) {
                cout << "(" << i+1 << "-";
                inTM = 1;
                nbTM--;
            }
            else if ((inTM==1) && (TM[i+1]==0)) {
                cout << i+1 << ") ";
                inTM = 0;
            }
        }
        if (nbTM) cout << "ERROR";
        cout << endl;
    }
    else { // outputFormat == 'L'
        cout << endl;
        for (int i = 0 ; i < nbAA ; i++ ) {
            cout << setw(4) << i+1 << "  ";
            cout << seqData[i] << "  ";
            cout << TM[i] << "  ";
            if (includeSPInfo) cout << TMInSwiss[i];
            cout << endl;
        }
    }        
}

bool ParseArg(int argc, char *argv[]) {
// parsing the arguments of the command line

    if ( argc < 2 ) {
      cerr << "PRED-TMR: at least one file name must be specified" << endl;
      return FALSE;
    }

    for (int i = 1 ; i < argc ; i++ ) {
        switch (argv[i][0] ) {
            case '-': {
                switch (argv[i][1]) {
                    case 'C' :
                    case 'c' : {
                        outputFormat = 'C';
                        break;
                    }
                    case 'R' :
                    case 'r' : {
                        outputFormat = 'R';
                        break;
                    }
                    case 'L' :
                    case 'l' : {
                        outputFormat = 'L';
                        break;
                    }
                    case 'T' :
                    case 't' : {
                        outputFormat = 'T';
                        break;
                    }
                    case 'O' :
                    case 'o' : {
                        showObservedTM = TRUE;
                        break;
                    }
                    case 'H' :
                    case 'h' : {
                        cerr << "PRED-TMR: usage" << endl;
                        cerr << "PRED-TMR [-[C|L|R]] [-O] [-H] <file names...>" << endl;
                        cerr << "-C = display result in condensed form (default)" << endl;
                        cerr << "-L = display result as a list" << endl;
                        cerr << "-R = display result by row" << endl;
                        cerr << "-T = display result in table mode" << endl;
                        cerr << "-O = show observed transmembrane regions" << endl;
                        cerr << "-H = print usage (this message)" << endl;
                        return FALSE;
                    }
                    default : {
                        cerr << "PRED-TMR: unknown option '" << argv[i][1] << "'" << endl;
                        return FALSE;
                    }

                }
                break;
            }
            default: {
                fileName[nbInputFile++] = i;
                break;
            }
        }
    }
    return 1;
}

bool ParseArg() {
    includeSPInfo = FALSE;
    char* seqText;
    cgidata cgi;
    for (int i=0; i< cgi.count(); i++) {

        if (strncmp(cgi.name(i), "Call_CoPreTHi",13) == 0) {
            call_CoPreTHi = TRUE;
        }
        if (strncmp(cgi.name(i), "Reset",5) == 0) {
            clearForm = TRUE;
        }
        else if (strcmp(cgi.name(i), "seq_name") == 0) {
            if (strlen(cgi.value(i)) > 0 ) strcpy(seqName, cgi.value(i));
        }
        else if (strcmp(cgi.name(i), "seq_data") == 0) {
			if (strlen(cgi.value(i)) > 10000) {
				cout << "Content-type: text/html" << endl << endl;
				cout << "Sequence must be less than 10000 character long !" << endl;
				return FALSE;
			}
            seqText = new char[strlen(cgi.value(i))+1];
            strcpy(seqText, cgi.value(i));
        }
        else if (strcmp(cgi.name(i), "show_obsTM") == 0) {
            if (strcmp(cgi.value(i), "yes") == 0 ) {
                showObservedTM = TRUE;
            }
        }
    }
    if (clearForm) {
        cout << "Content-type: text/html" << endl << endl;
        system("cat PRED-TMR.gab");
        return FALSE;
    }
    // write the sequence in a file called 'tempfile'
    ofstream tmpfile;
    do {
        tmpfile.clear();
        tmpfile.open("tempfile", ios::out | ios::noreplace);
    }
    while (tmpfile.fail());
    tmpfile << seqText;
    tmpfile << endl;
    tmpfile.close();

    delete[] seqText;
    return TRUE;
}

bool rating(float& Qa, float& Ca, float& SM) {
    float goodPred = 0;     // w
    float goodNonPred = 0;  // x
    float missedPred = 0;   // y
    float wrongPred = 0;    // z

    int inPredTM = 0;
    int inObsTM = 0;

    // to compute the correlation factor M & C
    // from the paper :
    // "Prediction of transmembrane a-helices in prokaryotic membrane proteins:
    // the dense alignment surface method"
    int nbObsTM = 0;  // number of experimental TM segments
    int nbExpM = 0;  // number of experimental TM segment that overlap with prediction
    int nbPredTM = 0; // number of predicted TM segments
    int nbPredM = 0; // number of predicted TM segment that overlap with experimental

    int expMatch = 0;
    int predMatch = 0;
    int total = 0;
    for (int i = 0 ; i < nbAA ; i++ ) {
        if ((TM[i] == 1) && (TMInSwiss[i] == 1)) goodPred++;
        if ((TM[i] == 0) && (TMInSwiss[i] == 0)) goodNonPred++;
        if ((TM[i] == 0) && (TMInSwiss[i] == 1)) missedPred++;
        if ((TM[i] == 1) && (TMInSwiss[i] == 0)) wrongPred++;

        if ((inObsTM == 0) && (TMInSwiss[i] == 1) ) {
            nbObsTM++;
            inObsTM = 1;
            expMatch = 0;
        }
        if ((inPredTM == 0) && (TM[i] == 1) ) {
            nbPredTM++;
            inPredTM = 1;
            predMatch = 0;
        }
        if ((inObsTM == 1) && (expMatch == 0) && (TM[i] == 1)) {
            expMatch = 1;
            nbExpM++;
        }

        if ((inPredTM == 1) && (predMatch == 0) && (TMInSwiss[i] == 1)) {
            predMatch = 1;
            nbPredM++;
        }

        if ((inObsTM == 1) && (TMInSwiss[i] == 0) ) inObsTM = 0;
        if ((inPredTM == 1) && (TM[i] == 0) ) inPredTM = 0;

        total++;
    }
    if (nbObsTM == 0) return FALSE;

    float pa = 100 * goodPred / (goodPred + missedPred);
    float pna = 100 * goodNonPred / (goodNonPred + wrongPred);
    Qa = (pa + pna)/2;
    float S = (goodPred+missedPred)/total;
    float P = (goodPred+wrongPred)/total;
    if ((P==0) || (S==0)) Ca = 0;
    else Ca = ((goodPred / total) - (P * S)) / fsqrt(P * S * (1 - P) * ( 1 - S ));

    float M = (nbObsTM == 0) ? 0 : (float)nbExpM/nbObsTM;
    float C = (nbPredTM == 0) ? 0 : (float)nbPredM/nbPredTM;
    SM = fsqrt(M*C);
    return TRUE;
}

void Generate_CoPreTHi_Page() {
    system("cat CoPreTHi_top.gab");
    cout << "<PARAM NAME=SEQ-NAME VALUE=";
    if (strlen(seqName) > 0) cout << seqName;
    else cout << "\"SEQUENCE NAME\"";
    cout << ">" << endl;
    cout << "<PARAM NAME=SEQ-LENGTH VALUE=" << nbAA << ">" << endl;
    int nbTM = 0;
    for (int i = 0 ; i < nbAA ; i++ ) {
        if (TM[i]==1) { // beginning of predicted TM region
            for (int j = i; j < nbAA ; j++ ) {
                if (TM[j]==0) break;
            }
            cout << "<PARAM NAME=PRED-TMR" << ++nbTM << " VALUE=";
            cout << "\"Transmem " << setw(2) << nbTM << "     " << setw(4) << i+1 << "   " << setw(4) << j << "\">" << endl;
            i=j;
        }
    }
    cout << "<PARAM NAME=NBPRED-TMR VALUE=" << nbTM << ">" << endl;
    nbTM = 0;
    for (i = 0 ; i < nbAA ; i++ ) {
        if (TMInSwiss[i]==1) { // beginning of observed TM region
            for (int j = i; j < nbAA ; j++ ) {
                if (TMInSwiss[j]==0) break;
            }
            cout << "<PARAM NAME=SWISSPROT" << ++nbTM << " VALUE=";
            cout << "\"FT   TRANSMEM   " << setw(4) << i+1 << "   " << setw(4) << j << "\">" << endl;
            i=j;
        }
    }
    cout << "<PARAM NAME=NBSWISSPROT VALUE=" << nbTM << ">" << endl;

    system("cat CoPreTHi_bottom.gab");
}
