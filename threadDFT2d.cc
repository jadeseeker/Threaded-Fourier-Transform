// Threaded two-dimensional Discrete FFT transform
// Yash Shah
// ECE 6122 Project 2


#include <iostream>
#include <string>
#include <math.h>
#include <stdlib.h>

#include "Complex.h"
#include "InputImage.h"

// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.

using namespace std;

Complex* imgData;
Complex* W ;
int imgWidth, imgHeight, imgSize;
int N, nThreads, flag;

int pThreads, count;
pthread_mutex_t countMutex;
bool* localSense;
bool globalSense;




// Function to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the
// number of points in the 1D transform.
unsigned ReverseBits(unsigned v)
{ //  Provided to students

  //N is the img width
  //No need to change
  unsigned n = N; // Size of array (which is even 2 power k value)
  unsigned r = 0; // Return value
   
  for (--n; n > 0; n >>= 1)
    {
      r <<= 1;        // Shift return value
      r |= (v & 0x1); // Merge in next bit
      v >>= 1;        // Shift reversal value
    }
  return r;
}

void bitRevereseOrder(Complex* input)
{

  Complex* output = new Complex[imgSize];
  for(int n = 0; n < N; n++)
  {
    int bufferSize = n*N; 
    Complex* temp = new Complex[N];  

    for(int i = 0; i < N; i++)
    {
      temp[i] = input[i + bufferSize];
    }
    for(int i = bufferSize; i < bufferSize + N; i++)
    {
      output[i] = temp[ReverseBits(i)];
    }
    delete [] temp;
  }
  
  for(int i = 0; i < N*N; i++)
  {
    input[i] = output[i];
  }
}


void MyBarrier_Init()
{

  count = pThreads;

  pthread_mutex_init(&countMutex,0);

  localSense = new bool[pThreads];
  
  for(int i = 0; i < pThreads; i++)
  {
    localSense[i] = true;
  }

  globalSense = true;  
}

// Each thread calls MyBarrier after completing the row-wise DFT
void MyBarrier(int myId) 
{

  localSense[myId] = !localSense[myId];
  pthread_mutex_lock(&countMutex);
  int myCount = count;
  count--;
  pthread_mutex_unlock(&countMutex);
  if(myCount == 1)
  {
    count = pThreads;
    globalSense = localSense[myId];
  }
  else
  {
    while(globalSense != localSense[myId]) { } // Spin
  }
    
}
                    
void Transform1D(Complex* h, int N)
{
  // Implement the efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assume even power of 2)
  Complex temp;
  for(int np = 2; np <= N; np = np*2)
  {
    for(int i = 0; i < N; i = i + np)
    {
      for(int j = 0; j < np/2; j++)
      {
        int offset = np/2;
        temp = h[i + j];
        h[i + j] = h[i + j] + W[j*N/np] * h[i + j + offset];
        h[i + j + offset] = temp - W[j*N/np] * h[i + j + offset];
      }
    }
  }

  if(flag == 0)
  {
    for(int i = 0; i < N; i++)
    {
      h[i].real = h[i].real/N;
      h[i].imag = h[i].imag/N;
    }
  }
}

void* Transform2DTHread(void* v)
{ // This is the thread startign point.  "v" is the thread number
  // Calculate 1d DFT for assigned rows
  // wait for all to complete
  // Calculate 1d DFT for assigned columns
  // Decrement active count and signal main if all complete

  // Assigning thread number
  unsigned long myId = (unsigned long)v; 

  // Performing 1D Transform
  int startRow = myId * N/nThreads;
  for(int n = 0; n < N/nThreads; n++)
  {
    int bufferSize = N*(startRow + n);  
    Transform1D(imgData + bufferSize, N);
  }

  // This thread has done its job; decrement the active count and see if all have finished
  MyBarrier(myId);
  
  return 0;
}

void calculateWeight(){
  //forward Transform
  if(flag == 1)
  {
    for(int i = 0; i < N/2; i++)
    {
      W[i] = Complex(cos(2 * M_PI * i / N), -sin(2 * M_PI * i / N));
    }
  }

  //inverse transform
  if(flag == 0)
  {
    for(int i = 0; i < N/2; i++)
    {
      W[i] = Complex(cos(2 * M_PI * i / N), sin(2 * M_PI * i / N));
    }
  }
}


void spawnThreads(){
  
  MyBarrier_Init();

  for (int i = 0; i < nThreads; ++i)
  {
    // Create the thread
    pthread_t pt; 
    pthread_create(&pt, 0, Transform2DTHread, (void*)i);
  }

}

void transpose(Complex* input)
{
  Complex* temp = new Complex[N*N];
  int pos = 0;
  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < N; j++)
    {
      temp[pos] = input[i + j*N];
      pos++;
    }
  } 

  for(int i = 0; i < N*N; i++)
  {
    input[i] = temp[i];
  }    
}

void trans1D(){

  bitRevereseOrder(imgData);
  //calculateWeight();
  spawnThreads();
  // Wait for all threads complete
  MyBarrier(nThreads);

}

void trans2D(){

  transpose(imgData);
  
  trans1D();

  transpose(imgData);
}

void Transform2D(const char* inputFN) 
{ 
  // Do the 2D transform here.
  // Create the helper object for reading the image
  InputImage image(inputFN);  
  // Create the global pointer to the image array data
  imgWidth = image.GetWidth();
  imgHeight = image.GetHeight();
  imgData = image.GetImageData();
  imgSize = imgHeight * imgWidth;

    // 1 for forward transform    0 for Inverse  
  flag = 1;

  //-----------------------------------------------------------------------------
  //        1D Transform
  //------------------------------------------------------------------------------
  
  calculateWeight();
  trans1D();

  cout<<"Generating Image File MyAfter1D.txt"<<endl;
  image.SaveImageData("MyAfter1D.txt", imgData, imgWidth, imgHeight);

  //-----------------------------------------------------------------------------
  //        2D Transform
  //------------------------------------------------------------------------------

  trans2D();

  cout<<"Generating Image File MyAfter2D.txt"<<endl;
  image.SaveImageData("MyAfter2D.txt", imgData, imgWidth, imgHeight); 


  //-----------------------------------------------------------------------------
  //        Inverse 2D Transform
  //------------------------------------------------------------------------------

  flag = 0;

  //Inverse 1D
  calculateWeight();
  trans1D();

  //Inverse 2D
  trans2D();

  for(int i = 0; i < N*N; i++)
  {
    if(imgData[i].real < 0.9  && imgData[i].imag < 0.9)
    {
      imgData[i].real = 0;
      imgData[i].imag = 0;
    }
  }

  cout<<"Generating Image File MyAfterInverse.txt"<<endl;
  image.SaveImageData("MyAfterInverse.txt", imgData, imgWidth, imgHeight);


}

int main(int argc, char** argv)
{

  N = 1024;
  nThreads = 16;
  
  W = new Complex[N/2];
  string fn("Tower.txt"); // default file name

  if (argc > 1 ){
    fn = string(argv[1]);  // if name specified on cmd line
    if(argc > 2){
      if(N % atoi(argv[2]) == 0){
        nThreads = atoi(argv[2]);
        //cout << "num threads: " << nThreads << endl;
      }
      else{
        cout << "Invalid Thread count : should be a factor of 1024 (dimension of img)" << endl;
        exit(0);
      }
    }
  }

  pThreads = nThreads + 1;
  // MPI initialization here
  Transform2D(fn.c_str()); // Perform the transform.
}  
  

  
