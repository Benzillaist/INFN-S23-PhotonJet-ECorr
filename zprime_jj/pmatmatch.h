#ifndef PMATMATCH_H
#define PMATMATCH_H

#include "butil.h"

int** matchJets(int nj_R, int nj_T, TVector3 recoVec3[], TVector3 truthVec3[], double minWeight) {

  // energy res weight: [0,1]
  double eResWeight = 0.2;
  // delta R weight: [0,1]
  double deltaRWeight = 0.8;

  int matchArrR[nj_T];
  int matchArrT[nj_T];
  int matchCount = 0;

  int minAxisLen = min(nj_R, nj_T);

  Double_t matchArrW[nj_R][nj_T];

  // cout << "nj_R: " << nj_R << " nj_T: " << nj_T << endl;

  for(int i = 0; i < nj_R; i++) {
    for(int j = 0; j < nj_T; j++) {
      matchArrW[i][j] = 1 - abs(( abs(eResWeight * (recoVec3[j] - truthVec3[i]).Mag() / truthVec3[i].Mag())) + (deltaRWeight * sqrt( pow(recoVec3[j].Eta() - truthVec3[i].Eta(), 2) + pow(recoVec3[j].Phi() - truthVec3[i].Phi(), 2) ) / TMath::Pi() ));
      // cout << abs(eResWeight * ((recoVec3[j] - truthVec3[i]).Mag()) / truthVec3[i].Mag())  << " " << (deltaRWeight * sqrt( pow(recoVec3[j].Eta() - truthVec3[i].Eta(), 2) + pow(recoVec3[j].Phi() - truthVec3[i].Phi(), 2) ) / TMath::Pi() ) << " " << matchArrW[i][j] << endl;
    }
  }

  int banR[nj_T];
  int banT[nj_T];
  int banCount = 0;
  for(int l = 0; l < minAxisLen; l++) {
    for(int i = 0; i < nj_T; i++) {
      Double_t maxWR = 0;
      Double_t maxWT = 0;
      int maxWRIndex = -1;

      // cout << "Test2.2" << endl;

      // cout << "nj_R: " << nj_R << " nj_T: " << nj_T << endl;

      for(int j = 0; j < nj_R; j++) {
        if(!contains(banR, banCount, j)) {
          if(matchArrW[j][i] > maxWR) {
            maxWR = matchArrW[j][i];
            maxWRIndex = j;
          }
        }
      }
      
      if(maxWRIndex != -1) {
        for(int j = 0; j < nj_T; j++) {
          if(!contains(banT, banCount, j)) {
            if(matchArrW[maxWRIndex][j] > maxWT) {
              maxWT = matchArrW[maxWRIndex][j];
            }
          }
        }
      }
      
      if((maxWT == maxWR) && (maxWT >= minWeight)) {
	// cout << "maxWR: " << maxWR << " maxWT: " << maxWT << " R: " << maxWRIndex << " T: " << i << endl;
        matchArrR[matchCount] = maxWRIndex;
        matchArrT[matchCount] = i;
        matchCount++;
        banR[banCount] = maxWRIndex;
        banT[banCount] = i;
        banCount++;
      }
    }
    
    if(banCount == minAxisLen) {
      break;
    }
  }

  int** ret = build2DArrI(3, nj_T);
  int size[nj_T];
  for(int i = matchCount; i < nj_T; i++) {
    matchArrR[i] = 0;
    matchArrT[i] = 0;
    size[i] = 0;
  }

  if(nj_T > 0) {
    size[0] = matchCount;
  }

  ret[0] = matchArrR;
  ret[1] = matchArrT;
  ret[2] = size;
  
  return ret;
}

#endif
