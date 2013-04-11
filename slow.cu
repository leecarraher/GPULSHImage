/*
Author: lee Carraher
 This file will start as a basic parallel centric framework 
for cuda and data streaming, which will eventually evolve into
supporting the leech decoder functionality and hopefully be 
able to directly borrow some functions from the previous decoders.
*/
// leech decoder Frameworks

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cuda.h>
#define blocksize 48
//this is a conversion for pointer arith. as cuda doesnt allow casting arrays
#define get(i,j,k) (j+(2*i+k)*4)

__device__ __constant__ int H6CodeWordsBin[16]= {88767,7138404,12026322,14300937,
                        7912737,1379322,13526604,10735767,
                        10429272,15897987,2750517,4476654,
                        15123654,9138717,5250987,4041072};



//this is an attempt to keep as much data as possible in the registers, such as the entire hexacode table
inline __device__  void H6Bin(char a,char b, char c,char* ret)
{
    unsigned int row = H6CodeWordsBin[a*4+b];
    row = ((row&0xFC0000) >>18)*(c==0) + ((row&0x03F000) >>12) * (c==1) + ((row&0x000FC0) >>6) * (c==2) + (row&0x00003F) *(c==3) ;
    ret[0] = (row&0x30)>>4;
    ret[1] = (row&0x0c)>>2;
    ret[2] =row&0x03;   
}
__device__ __constant__ int H6CodeWordsRevBin[16] ={235431, 11881674, 14219793, 7217532, 
                            5686002, 14686623, 9285444, 3896361, 
                            11091213, 2090592, 7491771, 12880854, 
                            16541784, 4895541, 2557422, 9559683};

inline __device__  void H6BinRev(char a,char b, char c,char ret[3])
{

    unsigned int row = H6CodeWordsRevBin[a*4+b];
    
    row = ((row&0xFC0000) >>18)*(c==0) + ((row&0x03F000) >>12) * (c==1) + ((row&0x000FC0) >>6) * (c==2) + (row&0x00003F) *(c==3) ;
    ret[0] = (row&0x30)>>4;
    ret[1] = (row&0x0c)>>2;
    ret[2] =row&0x03; 
}


__device__ __constant__ float qamPts[16][2] =
//000, 110 , 001, 111 aeven
 {{1.0, 7.0},{5.0, 7.0},{5.0, 3.0},{1.0, 3.0},
//010 100 011 101 aodd
{3.0, 5.0},{3.0, 1.0},{7.0, 1.0},{7.0, 5.0},
//000, 110 , 001, 111 beven
 {3.0, 7.0},{7.0, 7.0},{7.0, 3.0},{3.0, 3.0},
//010 100 011 101 bodd
 {5.0, 5.0},{5.0, 1.0},{1.0, 1.0},{1.0, 5.0}};

#define DIST(x1,x2, y1,y2) (((x1-x2)*(x1-x2)) + ((y1-y2) * (y1-y2)))


__device__ void QAM(float *r,int part,float *dijs,float *dijks,char *kparities){
    /*
        this function returns all of the pertinant information from the decoder such as minimum distances, nearest coset leader quadrant, and alternative k-parity distances
*/

    

    char i = 0;
    float dist000,dist110,dist001,dist111;
    char pos;
    char d;
    int qsection = 4*part;
    #pragma unroll 12
    for(;i<12;i++){

        dist000 = DIST(r[i*2],qamPts[qsection][0],r[i*2+1],qamPts[qsection][1]);
        dist110 = DIST(r[i*2],qamPts[1+qsection][0],r[i*2+1],qamPts[1+qsection][1]);
        dist001 = DIST(r[i*2],qamPts[2+qsection][0],r[i*2+1],qamPts[2+qsection][1]);
        dist111 = DIST(r[i*2],qamPts[3+qsection][0],r[i*2+1],qamPts[3+qsection][1]);
        
        //experimentating with removing if statements and thus branching
        //pos is the poition in the distance vectors
        pos = i*4+part;
        d = dist000<dist001;
        dijs[pos]=dist000*d+dist001*(!d);
        dijks[pos]=dist001*d+dist000*(!d);
        kparities[pos] = (!d);
        pos = i*4+(-part+3); //this result is from system of eqns
                            // satisfying 0*y +x =3  , 1*y +x =2
        d = dist110<dist111;
        dijs[pos]=dist110*d+dist111*(!d);
        dijks[pos]=dist111*d+dist110*(!d);
        kparities[pos] = qsection;//(!d);

    }
    

}

__device__ void blockConf(unsigned int part,float *dijs,float *mus,char *prefRep){
    /*
        computes the Z2 block confidence of the concatonated points projections onto GF4 characters
        this method combines previous versions branch removal and extends the 4 decoder parallel design
        to independently computer odd and even block confidences
    */



    //each two symbols is taken as a single character in GF4
    char d,pos,i=0;
    part = part%2;
    float s,t;
    #pragma unroll 6
    for(; i<6;i++){
    
    /* these are the part selector eqns
        0,15   3,12   10,5   6,9
        8,7    4,11   2,13   1,14
        --------------------------
       +8 -8  +1,-1  -8,+8  -5,+5

    */ 
        //0000 1111
        s = dijs[get(i,0,0)]+dijs[get(i,0,1)];
        t = dijs[get(i,3,0)]+dijs[get(i,3,1)];


        d = s<t;
        pos = i*4+part*24;
        mus[pos] = s*d+t*(!d);
        prefRep[pos] = (0+8*part)*d+(15-8*part)*(!d);


        //0011 1100 0 3 3 0
        s = dijs[get(i,0,0)]+dijs[get(i,3,1)];
        t = dijs[get(i,3,0)]+dijs[get(i,0,1)];
        
        d = s<t;
        pos = i*4+part*24+1;
        mus[pos] = s*d+t*(!d);
        prefRep[pos] = (3+1*part)*d+(12-1*part)*(!d);



        //1010 0101
        s = dijs[get(i,2,0)]+dijs[get(i,2,1)];
        t = dijs[get(i,1,0)]+dijs[get(i,1,1)];


        d = s<t;
        pos = i*4+part*24+2;
        mus[pos] = s*d+t*(!d);
        prefRep[pos] = (10-8*part)*d+(5+8*part)*(!d);


        //0110 1001
        s = dijs[get(i,1,0)]+dijs[get(i,2,1)];
        t = dijs[get(i,2,0)]+dijs[get(i,1,1)];
        
        d = s<t;
        pos = i*4+part*24+3;
        mus[pos] = s*d+t*(!d);
        prefRep[pos] = (6-5*part)*d+(9+5*part)*(!d);   
    }

}


////algo 2 from 1996 amrani/be'ery
__device__ float minH6(unsigned int part,char y[6],float* mus){
    /*
        this is the minimization over the hexacode funtion using the 2nd algorithm of  amrani and be'ery ieee may '96
    */
    
    if(part!=0)mus = (float *)&mus[6];
    float charwts[6] = {0,0,0,0,0,0};
    char i = 0;
    char d;
    char pos;
    //locate least weight characters in GF4
    #pragma unroll 6
    for(;i<6;i++)
    {
        char leastChar = 0;
        float leastwt = mus[i*4+0];
        pos = i*4;
        d = mus[pos+1]<leastwt;
        leastwt = mus[pos+1]*d+leastwt*(!d);
        leastChar = 1*d+leastChar*(!d);
        
        d = mus[pos+2]<leastwt;
        leastwt = mus[pos+2]*d+leastwt*(!d);
        leastChar = 2*d+leastChar*(!d);

        d = mus[pos+3]<leastwt;
        leastwt = mus[pos+3]*d+leastwt*(!d);
        leastChar = 3*d+leastChar*(!d);
        
        y[i] = leastChar;
        charwts[i]=leastwt;
    } 
    
    
    
    
    //test if equal
    char s[3];
    H6Bin(y[0],y[1],y[2],s);
    
    //how to get rid of this?
    //bool eval = (s[0]==y[3] && s[1]==y[4] && s[2]==y[5]);
    //return mus[0][y[0]]+mus[1][y[1]]+mus[2][y[2]]+mus[3][y[3]]+mus[4][y[4]]+mus[5][y[5]];


    //locate least reliable of the characters
    float leastreliablewt = charwts[0];
    char leastreliablechar = 0;
    
    d = charwts[1]>leastreliablewt;
    leastreliablewt = charwts[1]*d+ leastreliablewt*(!d);
    leastreliablechar = 1*d+ leastreliablechar*(!d);
    
    d = charwts[2]>leastreliablewt;
    leastreliablewt = charwts[2]*d+ leastreliablewt*(!d);
    leastreliablechar = 2*d+ leastreliablechar*(!d);
    
    
    char temp[6];
    
    temp[0] = y[0];
    temp[1] = y[1];
    temp[2] = y[2];
    temp[3] = y[3];
    temp[4] = y[4];
    temp[5] = y[5];
    
    //each computation of m_dist<minCodeWt will fail
    //temp will maintain y's which will be copied back to y
    //in the last step of this computation
    float m_dist,minCodeWt=1000.0;//*(1-eval);
    i = 0;
    #pragma unroll 4
    for(;i<4;i++)
    {
        y[leastreliablechar] = i;
        H6Bin(y[0],y[1],y[2],s);
        m_dist = mus[0*4+y[0]]+mus[1*4+y[1]]+mus[2*4+y[2]]+ 
                 mus[3*4+s[0]]+
                 mus[4*4+s[1]]+
                 mus[5*4+s[2]];


        d = m_dist < minCodeWt;
        minCodeWt = m_dist * d + minCodeWt*(!d);
        temp[0] = y[0]*d+temp[0]*(!d);
        temp[1] = y[1]*d+temp[1]*(!d);
        temp[2] = y[2]*d+temp[2]*(!d);
        temp[3] = s[0]*d+temp[3]*(!d);
        temp[4] = s[1]*d+temp[4]*(!d);
        temp[5] = s[2]*d+temp[5]*(!d);

    }


    
    //note y[0:3] are lost, we may need to back up

    

    //y2
    //locate the least reliable symbol in each
    leastreliablewt = charwts[3];
    leastreliablechar = 3;

    d = charwts[4]>leastreliablewt;
    leastreliablewt = charwts[4]*d+ leastreliablewt*(!d);
    leastreliablechar = 4*d+ leastreliablechar*(!d);
    
    d = charwts[5]>leastreliablewt;
    leastreliablewt = charwts[5]*d+ leastreliablewt*(!d);
    leastreliablechar = 5*d+ leastreliablechar*(!d);

    
    if(charwts[4]>leastreliablewt){
        leastreliablewt = charwts[4];
        leastreliablechar = 4;
    }
    if(charwts[5]>leastreliablewt){
        leastreliablewt = charwts[5];
        leastreliablechar = 5;
    }

    #pragma unroll 4
    for(i=0;i<4;i++)
    {
        
        y[leastreliablechar] = i;
        H6BinRev(y[3],y[4],y[5],s);
        m_dist = mus[0*4+s[0]]+
                 mus[1*4+s[1]]+
                 mus[2*4+s[2]]+
                 mus[3*4+y[3]]+mus[4*4+y[4]]+mus[5*4+y[5]];

        d = m_dist < minCodeWt;
        minCodeWt = m_dist * d + minCodeWt*(!d);
        temp[0] = s[0]*d+temp[0]*(!d);
        temp[1] = s[1]*d+temp[1]*(!d);
        temp[2] = s[2]*d+temp[2]*(!d);
        temp[3] = y[3]*d+temp[3]*(!d);
        temp[4] = y[4]*d+temp[4]*(!d);
        temp[5] = y[5]*d+temp[5]*(!d);

    }    


    //requires a deep copy here
    y[0] = temp[0];
    y[1] = temp[1];
    y[2] = temp[2];
    y[3] = temp[3];
    y[4] = temp[4];
    y[5] = temp[5];

    
    //this will return 0.0 distance if eval is true
    return minCodeWt;//*(eval-1);
}

__device__ float hparity(float weight,char hexword[6],char *prefReps,float *dijs,char oddFlag,unsigned int *codeword){
    /*
        here we are resolving the h-parity. which requres that the overall least significant bit parities equal the 
        bit parities of each projected GF4 block. aka column parity must equal 1st row parity
    */
    char parity= 0;
    char i =0;
    //#pragma unroll 6
    for(i=0;i<6;i++){
        parity = parity + (prefReps[i*4+hexword[i]]>7);//this should be the highest order bit  
        *codeword = *codeword + (prefReps[i*4+hexword[i]]<<(i*4));//should scoot for each 4bits, ok
    }
    
    //unrolled myself = unaligned memory access at yE line 840 something
    /*
    parity = parity + (prefReps[hexword[0]]>7);
    *codeword = *codeword + (prefReps[hexword[0]]);
    parity = parity + (prefReps[4+hexword[1]]>7);
    *codeword = *codeword + (prefReps[4+hexword[1]]<<(4));
    parity = parity + (prefReps[8+hexword[2]]>7);
    *codeword = *codeword + (prefReps[8+hexword[2]]<<(8));
    parity = parity + (prefReps[12+hexword[3]]>7);
    *codeword = *codeword + (prefReps[12+hexword[3]]<<(12));
    parity = parity + (prefReps[16+hexword[4]]>7);
    *codeword = *codeword + (prefReps[16+hexword[4]]<<(16));
    parity = parity + (prefReps[20+hexword[5]]>7);
    *codeword = *codeword + (prefReps[20+hexword[5]]<<(20));
    */




    unsigned int temp = *codeword;
    float leastwt = 1000.0;
    unsigned int least = 0;
    
    float deltaX;
    char idx1,idx2,idxComp1,idxComp2,proj;
    i = 0;
    char d;
    //walk along the codeword again
    #pragma unroll 6
    for(;i<6;i++){
        //bitwise method for complementing coordinates and evaluating dijs positions        
        proj = temp&0xf;//grab lower order 4 bits, 1111
        idx2 = proj&0x3;//grab second set of 2 bits ,    0011
        idx1 = (proj&0xc)>>2;//grab first set of 2 bits, 1100
        idxComp1 =idx1^0x3;//complement bits ^11
        idxComp2 =idx2^0x3;//complement bits ^11
        deltaX = (dijs[get(i,idxComp1,0)] + dijs[get(i,idxComp2,1)]) - (dijs[get(i,idx1,0)] + dijs[get(i,idx2,1)]);


        d = deltaX < leastwt;
        leastwt = deltaX*d+leastwt*(!d);
        least = i*d+least*(!d);
        /*
        if (deltaX < leastwt){
            leastwt = deltaX;
            least = i;
        }*/
        temp = temp>>4;//shift temp 4 bits
    }
    
    
    //#update the codeword and complement the correct section
    //print(*codeword);
    

        //dont break the codeword
    d = parity&1==oddFlag;
    /*if(parity&1==oddFlag)
        return weight;*/

    
    // (parity&1) - oddFlag, evaluates to 0 if parity is correct
    // thus replacing the branch with an extra computation
    weight = weight + leastwt*(!d) ;//* ((parity&1) - oddFlag);
    //bad ass bitwise method to complement the least_th 4 bit set

    *codeword= (*codeword ^ ( 0xF<<(least<<2) ))*(!d) + *codeword * (d);
    return weight;
}

__device__ float kparity(float weight,unsigned int codeword,char Btype, float *dijs,float *dijks,char* kparities){
    /*
        this last parity check assures that all A or B points have even/odd parity
    */



    unsigned int parity = 0;
    unsigned int i;
    unsigned int temp = codeword;

    #pragma unroll 12
    for(i=0 ;i <12;i++)
    {
        parity= parity+ kparities[i*4+(temp&3)];//&3 is bitmasking with 11, giving the low order bits
        temp=temp>>2;
    }
    
    
    //could return here but cuda does not like branching
    char d = parity&1 ==Btype;
    /*if(parity&1 ==Btype ){
        return weight;
    }*/
    
    

    float least = 1000 ;//*( parity&1 - Btype);
    float dif;
    char dd;
    #pragma unroll 12
    for(i=0 ;i <12;i++)
    {
        dif = dijks[i*4+(temp&3)]-dijs[i*4+(temp&3)];
        
        dd=dif < least;
        least = dif*dd + least*(!dd);
        //if(dif < least)least = dif;
        temp=temp>>2;

    }
    return least*(!d)+weight;
    
}





__device__ void decodePar4(float* r,int part,float *dijs, float *  dijks, char *  kparities, float * mus , char *prefreps , unsigned int *codewords,float *codewordweights,unsigned int *decoding )
{
     //unsigned int sharedoffset = 0;//threadIdx.x;


     QAM(r,part,  (float *)&dijs[part/2] ,(float *)&dijks[part/2],(char *)&kparities[part/2]);
     __syncthreads();
     //debug converts 32x8bytes into a single 32bit int, other 16 are lost, this is just a test
     //no sense in creating a special struct
     //for(i=0;i<32;i++)tmp = (tmp<<1)+ (unsigned int)kparities[sharedoffset][i];
     //return tmp;
     
     blockConf(part,(float *)&dijs[part/2],(float *)&mus[part/2],(char *)&prefreps[part/2]);
     
     
     __syncthreads();
    
    char y[6] = {0,0,0,0,0,0};
    float weight = minH6(part,y,(float *)&mus[part/2]);
    
    //is so far returning valid hexacodes
    if(part ==3 ){
    
        /*int i=0;
        float sum = 0.0;
        for(;i<24;i++)sum = sum + dijs[i*2]+r[i*2+1];*/
        *decoding = (unsigned int)(y[0]+(y[1]<<2)+(y[2]<<4)+(y[3]<<6)+(y[4]<<8)+(y[5]<<10));
    }
    
    
    unsigned int codeword = 0;
    weight = hparity(weight,y,(char *)&prefreps[part/2],(float *)&dijs[part/2],part%2,&codeword);//byref
    
    //these are not valid lattice points, but they may not be till after checking kparity----eesh thats not true
    //if(part==0){
    //    *decoding = codeword;
    //}
    
    
    
    weight = kparity(weight,codeword,(unsigned int)(part>2),(float *)&dijs[part/2],(float *)&dijks[part/2],(char *)&kparities[part/2]);
    //weights need to get stored in some shared region
    //    if(part==0){
    //    *decoding = codeword;
    //}
    
    /*only the part 0 decoder operates here to do comparisons
    it has to pull the datas gathered by the other simultaneous decoders
    sync theoretically is needed but the lock step cuda exec model may 
    force synchronization 
    */
    codewordweights[part]=weight;
    //codewords[part]=codeword;
    __syncthreads();
    
    /*if(part==0){
        float tmpCodewordWeight = codewordweights[0];
        *decoding =  codewords[0];
        
        if (tmpCodewordWeight > codewordweights[1]){
            tmpCodewordWeight = codewordweights[1];
            *decoding  =  codewords[1];
        }
        if (tmpCodewordWeight > codewordweights[2]){
            tmpCodewordWeight = codewordweights[2];
            *decoding  =  codewords[2];
        }
        if (tmpCodewordWeight > codewordweights[3]){
            tmpCodewordWeight = codewordweights[3];
            *decoding  =  codewords[3];
        }
    }*/
    
    
    
}

__global__ void compute(float *data,unsigned int* ret,int N)
{


   //mem for a block
   //each decoder has 4 threads so there must be threads_perblock/threads_perDecoder memory allocations
 __shared__ float dijs[blocksize/4][12][4][2];// A/B
 __shared__ float dijks[blocksize/4][12][4][2];//""
 __shared__ char kparities[blocksize/4][12][4][2];//""
 __shared__ float mus[blocksize/4][12][4][2];
 __shared__ char prefreps[blocksize/4][12][4][2];
 __shared__ unsigned int codewords[blocksize/4][4];
 __shared__ float codewordweights[blocksize/4][4];
 
  
  int id = threadIdx.x + blockIdx.x * blockDim.x;
  int partofquaddecoder = id%4;
  int decoderID = id/4  ;//0.1.2.3...blocksize ...2*blocksize ... N
  float* r = (float *)&data[decoderID*24];//each id is a 4 thread subdecoder
  
  decodePar4(r,partofquaddecoder, (float *)&dijs[decoderID] , (float *)  &dijks[decoderID] , (char *)  &kparities[decoderID] , (float *) &mus[decoderID] , (char *) &prefreps[decoderID] , (unsigned int *) &codewords[decoderID],(float*) &codewordweights[decoderID],&ret[decoderID]);
  
  //ret[decoderID]= (unsigned int)1000*dijs[decoderID][0][0][0] ;
  //partofquaddecoder = id%4
  //decoder = id/4
  /*
    DREWSTUBBS:
    ok so we need to first split up the data into blocks, first we were using shared memory and small block size as larger
    blocks resulted in register and shared limits. this was not as fast as using global memory (likely thanking the scheduler)
    now we need to figure out if we want to have shared blocks again (required for 4 simultaneous threads), and maybe have
    max_allowable_blocks/4 sub blocks, all using the same shared, or using seperate shared. it seems reasonable that we should
    test both. so how to roganize code to allow for this. also non-shared and shared block methods seem to peak around 32 blocks
    so perhaps a 32 blocks / 4 per decoding will result in a reasonable first round test.
  */

  
  
  
  }


int main(void)
{
  float *data_h;           // pointers to host memory
  unsigned int *ret_h;
  unsigned int *ret_d_h;            //point to host memory for device result
  float *data_d;           // pointer to device memory
  unsigned  int *ret_d;
  
                      //r0  1010 1010 1111 0000 1010 0101
      float r[4][12][2] ={{{ 7.0 , 5.0 }, { 3.0 , 1.0 }, { 3.0 , 1.0 }, { 3.0 , 1.0 }, 
                      { 1.0 , 3.0 }, { 1.0 , 3.0 }, { 1.0 , 7.0 }, { 5.0 , 3.0 }, 
                      { 7.0 , 5.0 }, { 7.0 , 5.0 }, { 7.0 , 1.0 }, { 7.0 , 1.0 }},
                      //r1 1011 0111 0010 1110 0100 1000
                      { { 3.0 , 1.0 }, { 5.0 , 7.0 }, { 7.0 , 1.0 }, { 5.0 , 7.0 }, 
                          { 5.0 , 3.0 }, { 7.0 , 5.0 }, { 1.0 , 3.0 }, { 3.0 , 1.0 }, 
                          { 7.0 , 1.0 }, { 5.0 , 3.0 }, { 7.0 , 5.0 }, { 5.0 , 3.0 } },
                      //r2  1111 1001 0011 1010 0000 1001
                      { { 3.0 , 3.0 }, { 3.0 , 3.0 }, { 5.0 , 1.0 }, 
                       { 5.0 , 5.0 }, { 3.0 , 7.0 }, { 3.0 , 3.0 },
                       { 5.0 , 1.0 }, { 1.0 , 5.0 }, { 3.0 , 7.0 },
                       { 3.0 , 7.0 }, { 1.0 , 5.0 }, { 5.0 , 5.0 }},
                      //r3  1011 0100 0111 0111 0100 0100
                      { { 5.0 , 1.0 }, { 3.0 , 3.0 }, { 1.0 , 1.0 }, { 3.0 , 7.0 },
                       { 1.0 , 1.0 }, { 7.0 , 7.0 }, { 5.0 , 5.0 }, { 3.0 , 3.0 }, 
                       { 1.0 , 1.0 }, { 7.0 , 3.0 }, { 1.0 , 1.0 }, { 3.0 , 7.0 }}}; 


  
  
  int i, N = 1000;
  size_t sizeIn = N*sizeof(float)*24;
  size_t sizeOut = N*sizeof(unsigned int);
  
  // allocate arrays on host
  data_h = (float *)malloc(sizeIn);
  ret_h = (unsigned  int *)malloc(sizeOut);
  ret_d_h = (unsigned  int *)malloc(sizeOut);
  
  srand ( time(NULL) );
  // Initialize host array and copy it to CUDA device
  for (int i=0; i<N; i++) 
  {
  
    for(int j = 0;j<12;j++)
    {
    
        data_h[(i*24)+(j*2)] = r[i%4][j][0];
        data_h[(i*24)+(j*2)+1] = r[i%4][j][1];
        //data_h[(i*24)+(j*2)] = float(rand()%9);
        //data_h[(i*24)+(j*2)+1] = float(rand()%9);
    }
  }
  
  // allocate array on device 
  cudaMalloc((void **) &data_d, sizeIn);
  cudaMalloc((void **) &ret_d, sizeOut);
  
  
  // copy data from host to device
  cudaMemcpy(data_d, data_h, sizeof(float)*N*24, cudaMemcpyHostToDevice);

  // Block Allocation
  int nBlocks = N/blocksize + (N%blocksize == 0?0:1);
  
  
  //Call kernel 
  compute<<< nBlocks, blocksize >>> (data_d,ret_d, N);
  
  // Retrieve result from device and store in b_h
  cudaMemcpy(ret_d_h, ret_d, sizeof(int)*N, cudaMemcpyDeviceToHost);
  //int j=0;
  for(i=0; i<4; i++){
    
    
    printf("%d\n",ret_d_h[i]);
    //show hex
    int j = 10;
    printf("r%d:",i%4);
    for(;j>=0;j-=2)printf("%d",(ret_d_h[i]>>j)&(3));
    

    //printf("%x  : ", ret_d_h[i]);
    //for(j=0;j<24;j++)printf("%f,",data_h[i*24+j]);
    printf("\n");
  }
  // cleanup
  free(data_h); free(ret_h); free(ret_d_h); cudaFree(data_d); cudaFree(ret_d); 
  //not so friendly remote upload
}

