/*Author: Lee Carraher
#Institution: University of Cincinnati, Computer Science Dept.


# this is a nearest lattice point decoder based on the hexacode based decoder of 
#Amrani, Be'ery IEEE Trans. on Comm. '96, with initial construction 
#from  Amrani, Be'ery,Vardy, Sun,Tilborg IEEE Info Thry'94

# the goal is to rewrite this algorithm in efficient C for cuda
# and eventual use as a Hashing Function
# for use in a Cuda Parallel Locality Hash Based Clustering algorthm
# additional implementation may include MPI/Cuda, and 
#anonymous offline data clusering


#-------------QAM Stuff ----------------------
# use a curtailed QAM for all positive signals
#  4 A000 B000 A110 B110
#  3 B101 A010 B010 A101 
#  2 A111 B111 A001 B001 
#  1 B011 A100 B100 A011 
#  0   1    2   3    4
# still gets rotated    \ 4 / 
#		                1 \/ 3
#                         /\ 
#	                    / 2 \	
# leech decoder uses a rotated Z2 lattice, so to find leading cosets
# just find the nearest point in 64QAM, A,B ; odd, even| to the rotated
# input vector
# rotate using the standard 2d rotation transform
#                      [cos x -sin x ]
#                  R = [sin x  cos x ]    cos(pi/4) = sin(pi/4)=1/sqrt(2)
# for faster C implementation use these binary fp constants
# 1/sqrt(2) = cc3b7f669ea0e63f ieee fp little endian
#           = 3fe6a09e667f3bcc ieee fp big endian
#           = 0.7071067811865475244008
#
#v' = v * R
# integer lattice
#
#  4 A000 B000 A110 B110 | A000 B000 A110 B110
#  3 B101 A010 B010 A101 | B101 A010 B010 A101
#  2 A111 B111 A001 B001 | A111 B111 A001 B001
#  1 B011 A100 B100 A011 | B011 A100 B100 A011
#    --------------------|---------------------
# -1 A000 B000 A110 B110 | A000 B000 A110 B110
# -2 B101 A010 B010 A101 | B101 A010 B010 A101
# -3 A111 B111v A001 B001 | A111 B111 A001 B001
# -4 B011 A100 B100 A011 | B011 A100 B100 A011
#even pts {000,110,111,001}
#odd  pts {010,101,100,011}
*/

#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>
#include <time.h>
#include "leechDecoder.c"

//#define DEBUG
//#define Algo1

static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))


#define HANDLE_NULL( a ) {if (a == NULL) { \
                            printf( "Host memory failed in %s at line %d\n", \
                                    __FILE__, __LINE__ ); \
                            exit( EXIT_FAILURE );}}


#define blocksize 32
//distance macro, inline support is different on different versions of cuda
#define DIST(x1,x2, y1,y2) (((x1-x2)*(x1-x2)) + ((y1-y2) * (y1-y2)))

//this is a conversion for pointer arith. as cuda doesnt allow casting arrays
#define get(i,j,k) (j+(2*i+k)*4)



#ifdef DEBUG
    static void print(int ret){
        int i;
        for(i=0;i<6;i++)
        {    
            printf("%d",(ret&8)>>3);
            printf("%d",(ret&4)>>2);
            printf("%d",(ret&2)>>1);
            printf("%d ",(ret&1));
            ret=ret>>4;
        } printf("\n");
    }
#endif



__device__ __constant__ int H6CodeWordsBin[16]= {88767,7138404,12026322,14300937,
                        7912737,1379322,13526604,10735767,
                        10429272,15897987,2750517,4476654,
                        15123654,9138717,5250987,4041072};



//this is an attempt to keep as much data as possible in the registers, such as the entire hexacode table
inline __device__  void H6Bin(char a,char b, char c,char* ret)
{
    unsigned int row = H6CodeWordsBin[a*4+b];
    /*switch(c){
        case 0:
            row = (row&0xFC0000) >>18;
        break;
        case 1: 
            row = (row&0x03F000) >>12;
        break;
        case 2:
            row = (row&0x000FC0) >>6;
        break;
        default:
            row = (row&0x00003F);
        break;
    }*/

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
    /*switch(c){
        case 0:
            row = (row&0xFC0000) >>18;
        break;
        case 1: 
            row = (row&0x03F000) >>12;
        break;
        case 2:
            row = (row&0x000FC0) >>6;
        break;
        default:
            row = (row&0x00003F);
        break;
    }*/



    ret[0] = (row&0x30)>>4;
    ret[1] = (row&0x0c)>>2;
    ret[2] =row&0x03; 
}



//000, 110 , 001, 111
__device__ __constant__ float evenAPts[4][2] = {{1.0, 7.0},{5.0, 7.0},{5.0, 3.0},{1.0, 3.0}};
//010 100 011 101
__device__ __constant__ float oddAPts[4][2]  ={{3.0, 5.0},{3.0, 1.0},{7.0, 1.0},{7.0, 5.0}};
//000, 110 , 001, 111
__device__ __constant__ float evenBPts[4][2] = {{3.0, 7.0},{7.0, 7.0},{7.0, 3.0},{3.0, 3.0}};
//010 100 011 101
__device__ __constant__ float oddBPts[4][2]  = {{5.0, 5.0},{5.0, 1.0},{1.0, 1.0},{1.0, 5.0}};



__device__ void QAM(float *r, float evenPts[4][2],float oddPts[4][2],float *dijs,float *dijks,char *kparities){
    /*
        this function returns all of the pertinant information from the decoder such as minimum distances, nearest coset leader quadrant, and alternative k-parity distances
    

    #these maps are seperated into the quadrants of a cartesian plane
    #now we gotta order these properly
      

    #another simple fix is that the quadrants of QAM be abstractly defined, and the -,+ of order
    #pairs be used to tile the generalized 16bit qam, besides this has to be done anyway so we
    #can get out the real number coordinates in the end
     */

    //the closest even-type Z2 lattice point is used as the 
    //coset representatives for all points, not currently used
    //quadrant = [0 for k in range(12)]

//change to float4 special datatype
//supposed to have better memory colescing

    

    char i = 0;
    float dist000,dist110,dist001,dist111,dist010,dist100,dist011,dist101;
    char pos;
    char d;
    #pragma unroll 12
    for(;i<12;i++){

        dist000 = DIST(r[i*2],evenPts[0][0],r[i*2+1],evenPts[0][1]);
        dist110 = DIST(r[i*2],evenPts[1][0],r[i*2+1],evenPts[1][1]);
        dist001 = DIST(r[i*2],evenPts[2][0],r[i*2+1],evenPts[2][1]);
        dist111 = DIST(r[i*2],evenPts[3][0],r[i*2+1],evenPts[3][1]);
        
        //experimentating with removing if statements and thus branching
        pos = i*4;
        d = dist000<dist001;
        dijs[pos]=dist000*d+dist001*(!d);
        dijks[pos]=dist001*d+dist000*(!d);
        kparities[pos] = (!d);
        /*if(dist000<dist001)
        {
             dijs[i*4]=dist000;
             dijks[i*4]=dist001;
             kparities[i*4] = 0;
        }
        else{
             dijs[i*4]=dist001;
             dijks[i*4]=dist000;
             kparities[i*4]= 1;
        }*/
        pos = i*4+3;
        d = dist110<dist111;
        dijs[pos]=dist110*d+dist111*(!d);
        dijks[pos]=dist111*d+dist110*(!d);
        kparities[pos] = (!d);


        /*if(dist110<dist111){
             dijs[i*4+3]=dist110;
             dijks[i*4+3]=dist111;
             kparities[i*4+3] = 0;
        }
        else{
             dijs[i*4+3]=dist111;
             dijks[i*4+3]=dist110;
             kparities[i*4+3] = 1;
        }*/
        //quadrant[i] = 0


        //min over odds
        dist010 = DIST(r[i*2],oddPts[0][0],r[i*2+1],oddPts[0][1]);
        dist100 = DIST(r[i*2],oddPts[1][0],r[i*2+1],oddPts[1][1]);
        dist011 = DIST(r[i*2],oddPts[2][0],r[i*2+1],oddPts[2][1]);
        dist101 = DIST(r[i*2],oddPts[3][0],r[i*2+1],oddPts[3][1]);
        
        pos = i*4+1;
        d = dist010<dist011;
        dijs[pos]=dist010*d+dist011*(!d);
        dijks[pos]=dist011*d+dist010*(!d);
        kparities[pos] = (!d);

        
        /*if (dist010<dist011){
             dijs[i*4+1]=dist010;
             dijks[i*4+1]=dist011;
             kparities[i*4+1] = 0;
        }
        else{
             dijs[i*4+1]=dist011;
             dijks[i*4+1]=dist010;
             kparities[i*4+1] = 1;   
        }*/

        pos = i*4+2;
        d = dist100<dist101;
        dijs[pos]=dist100*d+dist101*(!d);
        dijks[pos]=dist101*d+dist100*(!d);
        kparities[pos] = (!d);
        
        /*if (dist100<dist101){
             dijs[i*4+2]=dist100;
             dijks[i*4+2]=dist101;
             kparities[i*4+2] = 0;
        }
        else{
             dijs[i*4+2]=dist101;
             dijks[i*4+2]=dist100;
             kparities[i*4+2] = 1;
        }*/

    }

}




__device__ void blockConf(float *dijs,float *muEs,float *muOs,char *prefRepE,char *prefRepO){
    /*
        computes the Z2 block confidence of the concatonated points projections onto GF4 characters
        note
    */

    //each two symbols is taken as a single character in GF4
    char d,pos,i=0;
    
    float s,t;
    #pragma unroll 6
    for(; i<6;i++){
        
        //0000 1111
        s = dijs[get(i,0,0)]+dijs[get(i,0,1)];
        t = dijs[get(i,3,0)]+dijs[get(i,3,1)];


        d = s<t;
        pos = i*4+0;
        muEs[pos] = s*d+t*(!d);
        prefRepE[pos] = 0*d+15*(!d);

        /*if(s<t){
            muEs[i*4+0] = s;
            prefRepE[i*4+0] = 0;//[0,0,0,0]
        }
        else{
            muEs[i*4+0] = t;
            prefRepE[i*4+0] = 15;//[1,1,1,1]
        }*/

        //0011 1100 0 3 3 0
        s = dijs[get(i,0,0)]+dijs[get(i,3,1)];
        t = dijs[get(i,3,0)]+dijs[get(i,0,1)];
        
        d = s<t;
        pos = i*4+1;
        muEs[pos] = s*d+t*(!d);
        prefRepE[pos] = 3*d+12*(!d);


        /*if(s<t){
            muEs[i*4+1] = s;
            prefRepE[i*4+1] = 3;//[0,0,1,1]
        }
        else{
            muEs[i*4+1] = t;
            prefRepE[i*4+1] = 12;//[1,1,0,0]
        }*/


        //1010 0101
        s = dijs[get(i,2,0)]+dijs[get(i,2,1)];
        t = dijs[get(i,1,0)]+dijs[get(i,1,1)];


        d = s<t;
        pos = i*4+2;
        muEs[pos] = s*d+t*(!d);
        prefRepE[pos] = 10*d+5*(!d);

        /*
        if (s<t){
            muEs[i*4+2] = s;
            prefRepE[i*4+2] = 10;//[1,0,1,0]
            }
        else{
            muEs[i*4+2] = t;
            prefRepE[i*4+2] = 5;//[0,1,0,1]
        }*/

        //0110 1001
        s = dijs[get(i,1,0)]+dijs[get(i,2,1)];
        t = dijs[get(i,2,0)]+dijs[get(i,1,1)];
        
        d = s<t;
        pos = i*4+3;
        muEs[pos] = s*d+t*(!d);
        prefRepE[pos] = 6*d+9*(!d);


        /*if(s<t){
            muEs[i*4+3] = s;
            prefRepE[i*4+3] =6;// [0,1,1,0]
        }
        else{
            muEs[i*4+3] = t;
            prefRepE[i*4+3] = 9;//[1,0,0,1]
        }*/


    //this operation could be parallel, but probably doesnt need to be

        //1000 0111
        s = dijs[get(i,2,0)]+dijs[get(i,0,1)];
        t = dijs[get(i,1,0)]+dijs[get(i,3,1)];

        d = s<t;
        pos = i*4+0;
        muOs[pos] = s*d+t*(!d);
        prefRepO[pos] = 8*d+7*(!d);

        /*
        if(s<t){
            muOs[i*4+0] = s;
            prefRepO[i*4+0] = 8;//[1,0,0,0]
        }
        else{
            muOs[i*4+0] = t;
            prefRepO[i*4+0] = 7;//[0,1,1,1]
        }*/

        //0100 1011
        s = dijs[get(i,1,0)]+dijs[get(i,0,1)];
        t = dijs[get(i,2,0)]+dijs[get(i,3,1)];
        
        d = s<t;
        pos = i*4+1;
        muOs[pos] = s*d+t*(!d);
        prefRepO[pos] = 4*d+11*(!d);

        /*
        if (s<t){
            muOs[i*4+1] = s;
            prefRepO[i*4+1] = 4;//[0,1,0,0]
        }
        else{
            muOs[i*4+1] = t;
            prefRepO[i*4+1] = 11;//[1,0,1,1]
        }*/

        //0010 1101
        s = dijs[get(i,0,0)]+dijs[get(i,2,1)];
        t = dijs[get(i,3,0)]+dijs[get(i,1,1)];


        d = s<t;
        pos = i*4+2;
        muOs[pos] = s*d+t*(!d);
        prefRepO[pos] = 2*d+13*(!d);

        /*
        if(s<t){
            muOs[i*4+2] = s;
            prefRepO[i*4+2] =2;// [0,0,1,0]
        }
        else{
            muOs[i*4+2] = t;
            prefRepO[i*4+2] = 13;//[1,1,0,1]
        }*/

        //0001 1110
        s = dijs[get(i,0,0)]+dijs[get(i,1,1)];
        t = dijs[get(i,3,0)]+dijs[get(i,2,1)];
        
        d = s<t;
        pos = i*4+3;
        muOs[pos] = s*d+t*(!d);
        prefRepO[pos] = 1*d+14*(!d);

        /*
        if(s<t){
            muOs[i*4+3] = s;
            prefRepO[i*4+3] = 1;//[0,0,0,1]
        }
        else{
            muOs[i*4+3] = t;
            prefRepO[i*4+3] = 14;//[1,1,1,0]
        }*/
    }

}




////algo 2 from 1996 amrani/be'ery
__device__ float minH6(char y[6],float* mus){
    /*
        this is the minimization over the hexacode funtion using the 2nd algorithm of  amrani and be'ery ieee may '96
    */
    
    
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



        /*
        if(mus[i*4+1]<leastwt){
            leastwt = mus[i*4+1];
            leastChar = 1;
        }

        if(mus[i*4+2]<leastwt){
            leastwt = mus[i*4+2];
            leastChar = 2;
        }

        if(mus[i*4+3]<leastwt){
            leastwt = mus[i*4+3];
            leastChar = 3;
        }
        */
        
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

    /*
    if(charwts[1]>leastreliablewt){
        leastreliablewt = charwts[1];
        leastreliablechar = 1;
    }
    if(charwts[2]>leastreliablewt){
        leastreliablewt = charwts[2];
        leastreliablechar = 2;
    }*/

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

        /*if(m_dist < minCodeWt){
            minCodeWt = m_dist;
            temp[0] = y[0];temp[1] = y[1];temp[2] = y[2];
            temp[3] = s[0];
            temp[4] = s[1];
            temp[5] = s[2];
        }*/
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

        /*
        if(m_dist < minCodeWt){
            minCodeWt = m_dist;
            
            temp[0] = s[0];
            temp[1] = s[1];
            temp[2] = s[2];           
            temp[3] = y[3];temp[4] = y[4];temp[5] = y[5];
        }*/
    }    

    //printf("%d%d%d%d%d%d = %f\n",temp[0],temp[1],temp[2],temp[3],temp[4],temp[5],minCodeWt);
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



__device__ unsigned int decode(float* data)
{

     unsigned int sharedoffset = 0;//threadIdx.x;
     
     //__shared__ float dijs[blocksize][12][4] ;
     float dijs[1][12][4] ;
     //__shared__ float fdata[blocksize][12][4] ;
     float fdata[1][12][4] ;
     //__shared__ char cdata[blocksize][12][4];
     char cdata[1][12][4];
     float f2data[1][12][4] ;//hack fix to no redo QAM if we dont have to
     char c2data[1][12][4];

    QAM(data,evenAPts,oddAPts, (float *)&dijs[sharedoffset] ,(float *)&fdata[sharedoffset],(char *)&cdata[sharedoffset]);
    


    blockConf((float *)&dijs[sharedoffset],(float *)&f2data[0],(float *)&f2data[0][6],(char *)&c2data[0],(char *)&c2data[0][6]);
    
    

    // #####################Construct Hexacode Word ###################
    char yE[6] = {0,0,0,0,0,0};
    //float charwts[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    //constructHexWord(muEs,y,charwts);
    
    // #####################Minimize over the Hexacode ###################
    //float weight = minH6(y,muEs);
    float weightE = minH6(yE,(float *)&f2data[0]);
    

    


    #ifdef DEBUG
    printf("%d,%d,%d,%d,%d,%d = %f\n\n",yE[0],yE[1],yE[2],yE[3],yE[4],yE[5],weight);
    #endif
    
    //****chars = y = hexword ***** 
    unsigned int codewordE = 0;
    weightE = hparity(weightE,yE,(char *)&c2data[0],(float *)&dijs[sharedoffset],0,&codewordE);//byref
    
    
    
    //----------------A Odd Quarter Lattice Decoder----------------
    //constructHexWord(muOs,y,charwts);


    char yO[6] = {0,0,0,0,0,0};
    float weightO = minH6(yO,(float *)&f2data[0][6]);
    
    #ifdef DEBUG
    printf("%d,%d,%d,%d,%d,%d = %f\n\n",y[0],y[1],y[2],y[3],y[4],y[5],weight);
    #endif
    
    
    unsigned int codewordO = 0;
    weightO = hparity(weightO,yO,(char *)&c2data[0][6],(float *)&dijs[sharedoffset],1,&codewordO);//byref

    
    
    //gotta recompute dijk, and kparities
    //QAM(data,evenAPts,oddAPts, (float *)&dijs[sharedoffset] ,(float *)&fdata[sharedoffset],(char *)&cdata[sharedoffset]);
    
    
    weightE =kparity(weightE,codewordE,0,(float *)&dijs[sharedoffset],(float *)&fdata[sharedoffset],(char *)&cdata[sharedoffset]);
        


    
    weightO =kparity(weightO,codewordO,0, (float *)&dijs[sharedoffset], (float *)&fdata[sharedoffset],(char *)&cdata[sharedoffset]);

    float leastweight = weightE;
    unsigned int leastCodeword = codewordE;

    char d = weightO<leastweight;
    leastweight = weightO*d + leastweight*(!d);
    leastCodeword = codewordO*d + leastCodeword*(!d);
    /*if(weightO<leastweight)
    {
        leastweight = weightO;
        leastCodeword = codewordO;
    }*/

    //----------------H_24 Half Lattice Decoder for B points----------------
    
    
    QAM(data,evenBPts,oddBPts, (float *)&dijs[sharedoffset] ,(float *)&fdata[sharedoffset],(char *)&cdata[sharedoffset]);
    //QAM(r,evenBPts,oddBPts,dijs,dijks,kparities);
    
    
        blockConf((float *)&dijs[sharedoffset],(float *)&f2data[0],(float *)&f2data[0][6],(char *)&c2data[0],(char *)&c2data[0][6]);
    //blockConf(dijs,muEs,muOs,prefRepE,prefRepO);
    

    //----------------B Even Quarter Lattice Decoder----------------
    //constructHexWord(muEs,y,charwts);
    
    
    weightE = minH6(yE,(float *)&f2data[0]);
    
    
    #ifdef DEBUG
    printf("%d,%d,%d,%d,%d,%d = %f\n",y[0],y[1],y[2],y[3],y[4],y[5],weight);
    #endif
    codewordE = 0;
    weightE = hparity(weightE,yE,(char *)&c2data[0],(float *)&dijs[sharedoffset],0,&codewordE);//byref
    //weight = hparity(weight,y,prefRepE,dijs,0,&codeword);//byref

    

    //----------------B Odd Quarter Lattice Decoder----------------
    //constructHexWord(muOs,y,charwts);
    weightO = minH6(yO,(float *)&f2data[0][6]);
    #ifdef DEBUG
    printf("%d,%d,%d,%d,%d,%d = %f\n",y[0],y[1],y[2],y[3],y[4],y[5],weight);
    #endif

    codewordO = 0;
    weightO = hparity(weightO,yO,(char *)&c2data[0][6],(float *)&dijs[sharedoffset],1,&codewordO);//byref
      

    //QAM(data,evenBPts,oddBPts, (float *)&dijs[sharedoffset] ,(float *)&fdata[sharedoffset],(char *)&cdata[sharedoffset]);
        
    weightE =kparity(weightE,codewordE,1,(float *)&dijs[sharedoffset],(float *)&fdata[sharedoffset],(char *)&cdata[sharedoffset]);
    

    weightO =kparity(weightO,codewordO,1,(float *)&dijs[sharedoffset], (float *)&fdata[sharedoffset], (char *)&cdata[sharedoffset]);

    

    d = weightE<leastweight;
    leastweight = weightE*d + leastweight*(!d);
    leastCodeword = codewordE*d + leastCodeword*(!d);
    
    d = weightO<leastweight;
    leastweight = weightO*d + leastweight*(!d);
    leastCodeword = codewordO*d + leastCodeword*(!d);
    
    /*
    if(weightE<leastweight){
        leastweight = weightE;
        leastCodeword = codewordE;
    }
    
    
    if(weightO<leastweight){
        leastweight = weightO;
        leastCodeword = codewordO;
    }*/
    
    #ifdef DEBUG
    printf("least=%d weight = %f\n",leastCodeword,leastweight);
    #endif
    
    return leastCodeword;

}

__global__ void gpu_kernel(float* data ,unsigned int* ret, int N)
{
/*
        int srun =  blockDim.x*blocksize;

        int idx = threadIdx.x + blockIdx.x * blockDim.x;


        //while(idx<N){

            float* r = (float *)&data[idx*24];
            int fingerprint = decode(r);

            __syncthreads();

            ret[idx] = fingerprint;//fingerprint;
            idx = idx + srun;
*/

        //}


        int srun = N/blocksize + (N%blocksize?1:0);
        
        int idx = threadIdx.x + blockIdx.x * blockDim.x;

        //while(idx<srun){

            float* r = (float *)&data[(idx)*24];
            int fingerprint = decode(r);

            //__syncthreads();
            
            ret[idx] = fingerprint;//fingerprint;
            idx=idx+srun;

        
        //}
        
}

void seq_kernel(float* data ,unsigned int* ret, int N)
{
        int i = 0;

        float temp[12][2];
    
        for(;i<N;i++){
            float* r = (float *)&data[i*24];
            for(int j=0;j<12;j++){
                temp[j][0]=r[j*2];
                temp[j][1]=r[j*2+1];
            }
            ret[i] = decode(temp);
        }       
}



// main routine that executes on the host
int main(void)
{
  //must be less than max number of blocks 65535
  const int N = 65535*(blocksize);
  const int block_size = blocksize;
  //max number of blocks
  int n_blocks = 65535;//N/blocksize + (N%blocksize?1:0);//128;
  
  float *a_h, *a_d;  // Pointer to host & device arrays

  unsigned int *ret_h, *ret_d;// Pointer to return host & device arrays
  
  
  size_t size = 24*N * sizeof(float);
  a_h = (float *)malloc(size);        // Allocate array on host
  HANDLE_ERROR(cudaMalloc((void **) &a_d, size));   // Allocate array on device

  srand ( time(NULL) );
  // Initialize host array and copy it to CUDA device
  for (int i=0; i<N; i++) 
  {
    for(int j = 0;j<12;j++)
    {
        a_h[(i*24)+(j*2)] = float(rand()%9);
        a_h[(i*24)+(j*2)+1] = float(rand()%9);
    }
  }

  
  //allocate return data stream
  size_t sizeR = N * sizeof(unsigned int);
  unsigned int *retcmp_h;//compare fingerprints on cpu
  ret_h = (unsigned int *)malloc(sizeR);
  retcmp_h = (unsigned int *)malloc(sizeR);
  cudaMalloc((void **) &ret_d, sizeR); 

  HANDLE_ERROR(cudaMemcpy(a_d, a_h, size, cudaMemcpyHostToDevice));
  printf("Wrote Data to GPU: %i %i\n",n_blocks, block_size); 

  // Do calculation on device:
  printf("Calling GPU Kernel\n");
  int start = clock();
  gpu_kernel <<< n_blocks, block_size >>> (a_d,ret_d, N);

  // Retrieve result from device and store it in host array
  HANDLE_ERROR(cudaMemcpy(ret_h, ret_d, sizeof(unsigned int)*N, cudaMemcpyDeviceToHost));

  int gputime = clock()-start;

  printf("Calling Sequential Kernel\n");
  start = clock();

  seq_kernel(a_h,retcmp_h, N);

  int cputime = clock()-start;
  
  printf("Comparing Results\n");
  int count =0;
  for (int i=0; i<N;i++)
  { 
    if(ret_h[i]!=retcmp_h[i])count++;
	    //intf("\n%d = %d \n",ret_h[i],retcmp_h[i]);
	
  }
  printf("Incompatible Decodings: %d\n",count);
  printf("GPU Time:%d\n",gputime);
  printf("CPU Time:%d\n",cputime);
  printf("Speedup: %f\n",float(cputime)/float(gputime));
  // Cleanup
  free(a_h); cudaFree(a_d);
  free(ret_h); cudaFree(ret_d);
}



