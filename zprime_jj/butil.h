#ifndef BUTIL_H
#define BUTIL_H

// I'm too lazy to learn structs

bool charToBool(char i) {
  if(i == '0') {
    return false;
  }
  return true;
}

bool contains(int arr[], int l, int search) {
  for(int i = 0; i < l; i++) {
    if(arr[i] == search) {
      return true;
    }
  }
  return false;
}

bool contains(double arr[], int l, double search) {
  for(int i = 0; i < l; i++) {
    if(arr[i] == search) {
      return true;
    }
  }
  return false;
}

bool contains(TString arr[], int l, TString search) {
  for(int i = 0; i < l; i++) {
    if(arr[i] == search) {
      return true;
    }
  }
  return false;
}

int indexOf(int arr[], int l, int search) {
  for(int i = 0; i < l; i++) {
    if(arr[i] == search) {
      return i;
    }
  }
  return -1;
}

int indexOf(double arr[], int l, double search) {
  for(int i = 0; i < l; i++) {
    if(arr[i] == search) {
      return i;
    }
  }
  return -1;
}

int indexOf(TString arr[], int l, TString search) {
  for(int i = 0; i < l; i++) {
    if(arr[i] == search) {
      return i;
    }
  }
  return -1;
}

int** build2DArrI(int r, int c) {
  int** arr = new int*[r];
  for(int i = 0; i < r; i++) {
    arr[i] = new int[c];
  }
  return arr;
}

double** build2DArrD(int r, int c) {
  double** arr = new double*[r];
  for(int i = 0; i < r; i++) {
    arr[i] = new double[c];
  }
  return arr;
}

int* buildArrI(int l) {
  int* arr = new int[l];
  for(int i = 0; i < l; i++) {
    arr[i] = 0;
  }
  return arr;
}

double* buildArrD(int l) {
  double* arr = new double[l];
  for(int i = 0; i < l; i++) {
    arr[i] = 0;
  }
  return arr;
}

#endif
