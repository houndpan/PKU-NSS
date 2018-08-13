#include<gmp.h>
#include<stdio.h>
#include<time.h>
#include<string.h>
#include<math.h>

#define BOUND 1024      //需要多长的比特位数用于估计
#define CHANGE 512     //误差容忍上界，一旦预估误差到达便需要重新正规化向量

typedef struct Vector 
//定义向量组数据结构
{
    mpz_t f0;
    mpz_t f1;
    mpz_t g0;
    mpz_t g1;
    mpz_t F;
    mpz_t G;
    mpz_t f0_low;
    mpz_t f1_low;
    mpz_t g0_low;
    mpz_t g1_low;
    mpz_t F_low;
    mpz_t G_low;
    mpz_t f0_high;
    mpz_t f1_high;
    mpz_t g0_high;
    mpz_t g1_high;
    mpz_t temp[5];      //临时变量存储空间
    mpz_t cm[4];        //coeff matrix 是小f小g等量的系数矩阵
    mpz_t cF[4];        //coeff F   是大F对应的系数
    long powerofF;      //大F所借的2之幂次
    mpz_t cG[4];        //coeff G   是大G对应的系数
    long powerofG;      //大G所借的2之幂次
    long highlevel;     //高位估计变量所在的位置
    long lowlevel;      //低位估计变量的位置
}VECTOR;

void initVector(VECTOR*V)
//初始化向量组
{
    mpz_inits(V->F, V->f0, V->f0_high, V->f0_low, V->f1, V->f1_high, V->f1_low, V->F_low, V->G, V->g0, V->g0_high, V->g0_low, V->g1, V->g1_high, V->g1_low, V->G_low, NULL);
    mpz_inits(V->temp[0], V->temp[1], V->temp[2], V->temp[3], V->temp[4], NULL);
    mpz_inits(V->cm[0], V->cm[1], V->cm[2], V->cm[3], NULL);
    mpz_inits(V->cF[0], V->cF[1], V->cF[2], V->cF[3], NULL);
    mpz_inits(V->cG[0], V->cG[1], V->cG[2], V->cG[3], NULL);
    mpz_set_si(V->cm[0], 1);
    mpz_set_si(V->cm[1], 0);
    mpz_set_si(V->cm[2], 0);
    mpz_set_si(V->cm[3], 1);
    mpz_set_si(V->cF[0], 1);
    mpz_set_si(V->cF[1], 0);
    mpz_set_si(V->cF[2], 0);
    mpz_set_si(V->cF[3], 0);
    V->powerofF = 0;
    mpz_set_si(V->cG[0], 0);
    mpz_set_si(V->cG[1], 1);
    mpz_set_si(V->cG[2], 0);
    mpz_set_si(V->cG[3], 0);
    V->powerofG = 0;
    V->highlevel = 0;
    V->lowlevel = BOUND;
}

void clearVector(VECTOR*V)
//清除向量组占用的内存空间
{
    mpz_clears(V->F, V->f0, V->f0_high, V->f0_low, V->f1, V->f1_high, V->f1_low, V->F_low, V->G, V->g0, V->g0_high, V->g0_low, V->g1, V->g1_high, V->g1_low,  V->G_low, NULL);
    mpz_clears(V->temp[0], V->temp[1], V->temp[2], V->temp[3], V->temp[4], NULL);
    mpz_clears(V->cm[0], V->cm[1], V->cm[2], V->cm[3], NULL);
    mpz_clears(V->cF[0], V->cF[1], V->cF[2], V->cF[3], NULL);
    mpz_clears(V->cG[0], V->cG[1], V->cG[2], V->cG[3], NULL);
}

void Normal(VECTOR*V)
//对向量组正规化，详细计算出格点的具体值。
{
    long size = (long)(mpz_size(V->f0_high) * 64);
    if (size > 2 * BOUND) {V->highlevel = size - BOUND;}
    mpz_set(V->temp[4], V->F);
    mpz_mul(V->temp[0], V->F, V->cF[0]);    //计算F
    mpz_mul(V->temp[1], V->G, V->cF[1]);
    mpz_mul(V->temp[2], V->f1, V->cF[2]);
    mpz_mul(V->temp[3], V->g1, V->cF[3]);
    mpz_add(V->F, V->temp[0], V->temp[1]);
    mpz_add(V->temp[0], V->temp[2], V->temp[3]);
    mpz_add(V->F, V->F, V->temp[0]);
    mpz_tdiv_q_2exp(V->F, V->F, V->powerofF);

    mpz_mul(V->temp[0], V->temp[4], V->cG[0]);    //计算G
    mpz_mul(V->temp[1], V->G, V->cG[1]);
    mpz_mul(V->temp[2], V->f1, V->cG[2]);
    mpz_mul(V->temp[3], V->g1, V->cG[3]);
    mpz_add(V->G, V->temp[0], V->temp[1]);
    mpz_add(V->temp[0], V->temp[2], V->temp[3]);
    mpz_add(V->G, V->G, V->temp[0]);
    mpz_tdiv_q_2exp(V->G, V->G, V->powerofG);

    mpz_mul(V->temp[0], V->f0, V->cm[0]);    //计算f0与g0
    mpz_mul(V->temp[1], V->g0, V->cm[1]);
    mpz_mul(V->temp[2], V->f0, V->cm[2]);
    mpz_mul(V->temp[3], V->g0, V->cm[3]);
    mpz_add(V->f0, V->temp[0], V->temp[1]);
    mpz_add(V->g0, V->temp[2], V->temp[3]);

    mpz_mul(V->temp[0], V->f1, V->cm[0]);    //计算f1与g1
    mpz_mul(V->temp[1], V->g1, V->cm[1]);
    mpz_mul(V->temp[2], V->f1, V->cm[2]);
    mpz_mul(V->temp[3], V->g1, V->cm[3]);
    mpz_add(V->f1, V->temp[0], V->temp[1]);
    mpz_add(V->g1, V->temp[2], V->temp[3]);

    mpz_tdiv_r_2exp(V->F_low, V->F, V->lowlevel);   //计算F_low
    mpz_tdiv_r_2exp(V->G_low, V->G, V->lowlevel);   //计算G_low
    mpz_tdiv_r_2exp(V->f0_low, V->f0, V->lowlevel);   //计算f0_low
    mpz_tdiv_r_2exp(V->f1_low, V->f1, V->lowlevel);   //计算f1_low
    mpz_tdiv_r_2exp(V->g0_low, V->g0, V->lowlevel);   //计算g0_low
    mpz_tdiv_r_2exp(V->g1_low, V->g1, V->lowlevel);   //计算g1_low
    
    mpz_tdiv_q_2exp(V->f0_high, V->f0, V->highlevel);     //计算f0_high
    mpz_tdiv_q_2exp(V->f1_high, V->f1, V->highlevel);     //计算f1_high
    mpz_tdiv_q_2exp(V->g0_high, V->g0, V->highlevel);     //计算g0_high
    mpz_tdiv_q_2exp(V->g1_high, V->g1, V->highlevel);     //计算g1_high

    mpz_set_si(V->cm[0], 1);       //重新初始化系数矩阵
    mpz_set_si(V->cm[1], 0);
    mpz_set_si(V->cm[2], 0);
    mpz_set_si(V->cm[3], 1);
    mpz_set_si(V->cF[0], 1);
    mpz_set_si(V->cF[1], 0);
    mpz_set_si(V->cF[2], 0);
    mpz_set_si(V->cF[3], 0);
    V->powerofF = 0;
    mpz_set_si(V->cG[0], 0);
    mpz_set_si(V->cG[1], 1);
    mpz_set_si(V->cG[2], 0);
    mpz_set_si(V->cG[3], 0);
    V->powerofG = 0;
}


void Max_mpz(mpz_t rop, mpz_t op1, mpz_t op2)
//将op1与op2中比较大的数赋值给rop
{
    if (mpz_cmp(op1, op2) > 0) {       //rop = MAX{op1, op2}
        mpz_set(rop, op1);
    } else {
        mpz_set(rop, op2);
    };
}

int Max_cmp_pair(mpz_t op00, mpz_t op01, mpz_t op10, mpz_t op11, mpz_t max0, mpz_t max1, mpz_t temp0, mpz_t temp1) 
//比较两组数绝对值化后，最大数在哪个数组
//如果|op00|或|op01|最大，则返回1, 否则返回0
{
    mpz_abs(temp0, op00);   //temp0 = |op00|
    mpz_abs(temp1, op01);   //temp1 = |op01|
    Max_mpz(max0, temp0, temp1);    //max0 = MAX{|op00|, |op01|} 
    mpz_abs(temp0, op10);
    mpz_abs(temp1, op11);
    Max_mpz(max1, temp0, temp1);    //max1 = MAX{|op10|, |op11|}
    if (mpz_cmp(max0, max1) > 0){
        return 1;                   //if max0 > max1 then return 1; else return 0; fi;
    } else {
        return 0;
    }
}

int judgeone(mpz_t op00, mpz_t op01, mpz_t op10, mpz_t op11, mpz_t temp00, mpz_t temp01, mpz_t temp10, mpz_t temp11, mpz_t temp_f, mpz_t temp_g, mpz_t temp1, mpz_t temp2)
//判断d=1或d=-1是否可以最小化
{
    mpz_add(temp00, op00, op10);    //temp00 = f0 + g0
    mpz_add(temp01, op01, op11);    //temp01 = f1 + g1
    mpz_sub(temp10, op00, op10);    //temp10 = f0 - g0
    mpz_sub(temp11, op01, op11);    //temp11 = f1 - g1
    if (Max_cmp_pair(temp00, temp01, temp10, temp11, temp_f, temp_g, temp1, temp2) == 0){        //if phi(f+g) < phi(f-g)
        mpz_addmul_ui(temp10, op10, 4);     //temp10 = f0 + 3 * g0
        mpz_addmul_ui(temp11, op11, 4);     //temp11 = f1 + 3 * g1
        if (Max_cmp_pair(temp00, temp01, temp10, temp11, temp_f, temp_g, temp1, temp2) == 0) {return 1;}    //if phi(f+g) < phi(f+3g)
    } else {
        mpz_submul_ui(temp00, op10, 4);     //temp00 = f0 - 3 * g0
        mpz_submul_ui(temp01, op11, 4);     //temp01 = f1 - 3 * g1
        if (Max_cmp_pair(temp00, temp01, temp10, temp11, temp_f, temp_g, temp1, temp2) == 1) {return -1;}   //if phi(f-g) < phi(f-3g)
    }
    return 0;
}

int Minimize(mpz_t op00, mpz_t op01, mpz_t op10, mpz_t op11, mpz_t rop, mpz_t temp_f, mpz_t temp_g, mpz_t temp00, mpz_t temp01, mpz_t temp10, mpz_t temp11, mpz_t temp1,  mpz_t temp2)
//计算使得两组数组合后绝对值最小的奇数//Caulate odd number d which Minimum MAX{|op00 + d * op10|, |op01 + d*op11|}
//算法的原理是这样的，考虑两个函数方程y0(x) = |f0 + g0 * x| 和 y1(x) = |f1 + g1 * x|, 我们所求即MAX{y0, y1}在x取值为奇数时的最小值。也就是函数MAX{y0(x), y1(x)}的最小值。
//注意到内层两个函数的形状都是先减后增的线性函数，所以他们的最大值函数也一定是一个先减后增的函数，并且我们可以计算出x取何值时整体取最小值。
//当g0 和 g1同号时，最小值点处于方程f0 + g0 * x = -(f1 + g1 * x)的解处， 当g0和g1异号时，最小值点位于f0 + g0 * x = f1 + g1 * x的解处
//前者的解为-(f0 + f1)/(g0 + g1), 后者的解为(f1 - f0)/(g0 - g1)
//得到最小值点后我们再比较距离其最近的两个奇整数点位的大小，选取其值最小的d输出。
//除op**外，输入参数中的rop、temp_*等皆为便于计算使用的临时存储空间，算法上无实际意义
{  
    int d;
    if (mpz_sgn(op10) != -1) {
        if (mpz_sgn(op11) != -1){       //判断g1和g0是否同号
            mpz_add(temp_f, op00, op01); 
            mpz_neg(temp_f, temp_f);
            mpz_add(temp_g, op10, op11);
            mpz_fdiv_q(rop, temp_f, temp_g);    //rop是MAX{y0(x), y1(x)}在实数域上的最小值点向下取整的结果。
            d = mpz_get_si(rop);
        } else {            //其余选择分支同理
            mpz_sub(temp_f, op01, op00);
            mpz_sub(temp_g, op10, op11);
            mpz_fdiv_q(rop, temp_f, temp_g);
            d = mpz_get_si(rop);
        }
    } else {
        if (mpz_sgn(op11) != -1){
            mpz_sub(temp_f, op01, op00);
            mpz_sub(temp_g, op10, op11);
            mpz_fdiv_q(rop, temp_f, temp_g);
            d = mpz_get_si(rop);
        } else {
            mpz_add(temp_f, op00, op01);
            mpz_neg(temp_f, temp_f);
            mpz_add(temp_g, op10, op11);
            mpz_fdiv_q(rop, temp_f, temp_g);
            d = mpz_get_si(rop);
        }
    }
    if (d % 2 == 0) {
        d = d - 1;          //2|d时，d - 1 和 d + 1是距离最小值点最近的两个奇数。否则d 和 d + 2是距离最小值点最近的两个奇数。
    }
    mpz_set(temp00, op00);
    mpz_set(temp01, op01);
    mpz_mul_si(temp1, op10, d);
    mpz_add(temp00, temp00, temp1);         //temp00 = f0 + d * g0 
    mpz_mul_si(temp2, op11, d);
    mpz_add(temp01, temp01, temp2);         //temp01 = f1 + d * g1 
    mpz_set(temp10, temp00);
    mpz_set(temp11, temp01);
    mpz_mul_si(temp1, op10, 2);
    mpz_add(temp10, temp10, temp1);         //temp10 = f0 + (d + 2) * g0 
    mpz_mul_si(temp2, op11, 2);
    mpz_add(temp11, temp11, temp2);         //temp11 = f1 + (d + 2) * g1
    if (Max_cmp_pair(temp10, temp11, temp00, temp01, temp1, temp2, temp_f, temp_g) > 0){
        return d;
    } else {
        return d + 2;
    };
}

void Rational_Approximation_ex(char* Input, mpz_t output[2])
{
    long d;      //使得f与g组合后最小的奇数d
    long Length = strlen(Input);          //输入01串的长度
    char state = Input[0];               //初始化指针在输入字符串中的位置
    long i = 1;                          //标记指针的位置为i - 1
    int j;                                //循环变量下标
    double dv = 0;                      //预估误差位数
    long run = 2;                       //游程计数器
    VECTOR cycle;                       //向量组
    initVector(&cycle);             //初始化向量组

    //初始化计时器
    clock_t start, finish;
    double timegap;
    start = clock();
    
    //排除输入序列前缀的零
    while(state == '0'){        //读前若干个0，确定第一个非零数前有i-1个零
        i++;
        state = Input[i-1];
    };

    //临时存储变量动态内存的初始化
    mpz_t temp_0, temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7, temp_8;           //这些是程序中可能用到的一些中间变量的存储位置，避免频繁分配空间，加速程序
    mpz_inits(temp_0, temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7, temp_8, NULL);

    //根据前缀零的个数决定f，g，F，G的初始值
    mpz_set_si(cycle.f0,0);
    mpz_set_si(cycle.f1,2);
    mpz_ui_pow_ui(cycle.g0,2,i-1);
    mpz_set_si(cycle.g1,1);
    mpz_set_si(cycle.F,1);        //保证永远都有F = （f1 × alpha - f0）/ 2^{i} 为整数, 初始时有alpha = 2^{i-1}
    mpz_set_si(cycle.G,0);        //保证永远都有G = （g1 × alpha - g0）/ 2^{i} 为整数
    Normal(&cycle);
    //迭代有理逼近开始
    
    
    while(i < Length){

        //计时显示输出部分，查看当前进度
        if (i % 100000 == 0) {
            finish = clock();
            timegap = (double)(finish - start)/ CLOCKS_PER_SEC;
            printf("i = %ld, 目前耗时%f 秒\n", i, timegap);
        }

        //当a[i] = 1 时，alpha = alpha + 2^{i}
        //于是我们计算 F，G的更新：F = F + f1; G = G + g1;
        //此时 g1 × alpha - g0 = G * 2^{i}, 如果 g1 × alpha - g0 = 0 mod 2^{i+1} 即 G = 0 mod 2
        //于此计算出新的f和g之后修正F为 (F + d * G)/2 或（2F）/2, G同理。之后i++进入下一轮循环
        
        if (Input[i] == '1'){
            mpz_mul_2exp(temp_0, cycle.cm[0], cycle.powerofF);
            mpz_add(cycle.cF[2], cycle.cF[2], temp_0);
            //cycle.cF[2] = cycle.cF[2] + (int)pow(2, cycle.powerofF) * cycle.cm[0];
            mpz_mul_2exp(temp_0, cycle.cm[1], cycle.powerofF);
            mpz_add(cycle.cF[3], cycle.cF[3], temp_0);
            //cycle.cF[3] = cycle.cF[3] + (int)pow(2, cycle.powerofF) * cycle.cm[1];
            mpz_mul_2exp(temp_0, cycle.cm[2], cycle.powerofG);
            mpz_add(cycle.cG[2], cycle.cG[2], temp_0);
            //cycle.cG[2] = cycle.cG[2] + (int)pow(2, cycle.powerofG) * cycle.cm[2];
            mpz_mul_2exp(temp_0, cycle.cm[3], cycle.powerofG);
            mpz_add(cycle.cG[3], cycle.cG[3], temp_0);
            //cycle.cG[3] = cycle.cG[3] + (int)pow(2, cycle.powerofG) * cycle.cm[3];
            mpz_add(cycle.F_low, cycle.F_low, cycle.f1_low);
            mpz_add(cycle.G_low, cycle.G_low, cycle.g1_low);
        }

        //判断旧的g是否满足g1 × alpha - g0 = 0 mod 2^{i+1}
        //如果满足，那么g不变，f乘2, G除2, F不变
        //如果不满足，则依据f与g的值挑选较大者减去奇数倍的较小者来得出新的f、g以及F、G。具体算法细节参看论文。
        if ( mpz_divisible_2exp_p(cycle.G_low, 1) != 0 ) {

            mpz_mul_2exp(cycle.cm[0], cycle.cm[0], 1);
            //cycle.cm[0] = cycle.cm[0] * 2;
            mpz_mul_2exp(cycle.cm[1], cycle.cm[1], 1);
            //cycle.cm[1] = cycle.cm[1] * 2;
            cycle.powerofG++;

            mpz_mul_2exp(cycle.f0_low, cycle.f0_low, 1);
            mpz_mul_2exp(cycle.f1_low, cycle.f1_low, 1);
            mpz_mul_2exp(cycle.f0_high, cycle.f0_high, 1);
            mpz_mul_2exp(cycle.f1_high, cycle.f1_high, 1);
            mpz_tdiv_q_2exp(cycle.G_low, cycle.G_low, 1);

            dv++;
            run++;      //游程计数器加一
        } else if (Max_cmp_pair(cycle.f0_high, cycle.f1_high, cycle.g0_high, cycle.g1_high, temp_0, temp_1, temp_2, temp_3) != 0) {             //如果旧的g不满足，我们就要寻找新的g，这个分支是f在给定范数下大于g的情况。
            //在游程为零时（即上一比特刚刚纠错过），d大概率等于正负一，所以此时以judgeone方法判断d的具体值
            if (run < 1) {
                d = judgeone(cycle.f0_high, cycle.f1_high, cycle.g0_high, cycle.g1_high, temp_0, temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7);
                if (d == 0) d = Minimize(cycle.f0_high, cycle.f1_high, cycle.g0_high, cycle.g1_high, temp_0, temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7, temp_8);
            } else {
                d = Minimize(cycle.f0_high, cycle.f1_high, cycle.g0_high, cycle.g1_high, temp_0, temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7, temp_8);
            }
            //这是预估误差的增加值
            dv = dv + ceil(log((double)fabs(d)));
            
            mpz_mul_si(temp_0, cycle.cm[2], d);             //按照d修改系数矩阵
            mpz_add(temp_0, temp_0, cycle.cm[0]);
            //Temp_num[0] = cycle.cm[0] + d * cycle.cm[2];
            mpz_mul_si(temp_1, cycle.cm[3], d);
            mpz_add(temp_1, temp_1, cycle.cm[1]);
            //Temp_num[1] = cycle.cm[1] + d * cycle.cm[3];
            mpz_mul_2exp(cycle.cm[0], cycle.cm[2], 1);
            //cycle.cm[0] = 2 * cycle.cm[2];
            mpz_mul_2exp(cycle.cm[1], cycle.cm[3], 1);
            //cycle.cm[1] = 2 * cycle.cm[3];
            mpz_set(cycle.cm[2], temp_0);
            //cycle.cm[2] = Temp_num[0];
            mpz_set(cycle.cm[3], temp_1);
            //cycle.cm[3] = Temp_num[1];

            mpz_mul_si(temp_0, cycle.g0_low, d);            //按照d修改低位信息
            mpz_add(cycle.f0_low, cycle.f0_low, temp_0);
            mpz_mul_si(temp_1, cycle.g1_low, d);
            mpz_add(cycle.f1_low, cycle.f1_low, temp_1);
            mpz_mul_2exp(cycle.g0_low, cycle.g0_low, 1);
            mpz_mul_2exp(cycle.g1_low, cycle.g1_low, 1);
            mpz_swap(cycle.g0_low, cycle.f0_low);
            mpz_swap(cycle.g1_low, cycle.f1_low);

            mpz_mul_si(temp_0, cycle.g0_high, d);           //按照d修改高位信息
            mpz_add(cycle.f0_high, cycle.f0_high, temp_0);
            mpz_mul_si(temp_1, cycle.g1_high, d);
            mpz_add(cycle.f1_high, cycle.f1_high, temp_1);
            mpz_mul_2exp(cycle.g0_high, cycle.g0_high, 1);
            mpz_mul_2exp(cycle.g1_high, cycle.g1_high, 1);
            mpz_swap(cycle.g0_high, cycle.f0_high);
            mpz_swap(cycle.g1_high, cycle.f1_high);
            
            if (cycle.powerofG > cycle.powerofF) {           //按照d修改F与G的系数信息
                for (j = 0; j < 4; j++){
                    mpz_mul_2exp(cycle.cF[j], cycle.cF[j], cycle.powerofG - cycle.powerofF);
                    //cycle.cF[j] = cycle.cF[j] * (int)pow(2, cycle.powerofG - cycle.powerofF);
                }
                cycle.powerofF = cycle.powerofG;
            } else {
                for (j = 0; j < 4; j++){
                    mpz_mul_2exp(cycle.cG[j], cycle.cG[j], cycle.powerofF - cycle.powerofG);
                    //cycle.cG[j] = cycle.cG[j] * (int)pow(2, cycle.powerofF - cycle.powerofG);
                }
                cycle.powerofG = cycle.powerofF;
            }
            for(j = 0; j < 4; j++) {
                mpz_mul_si(temp_0, cycle.cG[j], d);
                mpz_add(cycle.cF[j], cycle.cF[j], temp_0);
                //cycle.cF[j] = cycle.cF[j] + d * cycle.cG[j];
            }
            for(j = 0; j < 4; j++){
                mpz_swap(cycle.cF[j], cycle.cG[j]);
                //Temp_num[0] = cycle.cF[j];
                //cycle.cF[j] = cycle.cG[j];
                //cycle.cG[j] = Temp_num[0];
            }
            cycle.powerofG = cycle.powerofG + 1;
            
            mpz_mul_si(temp_0, cycle.G_low, d);       //按照d修改F与G的低位信息
            mpz_add(cycle.F_low, cycle.F_low, temp_0);
            mpz_tdiv_q_2exp(cycle.F_low, cycle.F_low, 1);
            mpz_swap(cycle.F_low, cycle.G_low);

            run = 0;    //归零游程计数器
        } else {             //如果旧的g不满足，我们就要寻找新的g，这个分支是g在给定范数下大于f的情况。            
            if (run < 1) {
                d = judgeone(cycle.g0_high, cycle.g1_high, cycle.f0_high, cycle.f1_high, temp_0, temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7);
                if (d == 0) d = Minimize(cycle.g0_high, cycle.g1_high, cycle.f0_high, cycle.f1_high, temp_0, temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7, temp_8);
            } else {
                d = Minimize(cycle.g0_high, cycle.g1_high, cycle.f0_high, cycle.f1_high, temp_0, temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7, temp_8);
            }
            dv = dv + ceil(log((double)fabs(d)));
            
            mpz_mul_si(temp_0, cycle.cm[0], d);                       //按照d修改系数矩阵
            mpz_add(cycle.cm[2], cycle.cm[2], temp_0);
            //cycle.cm[2] = cycle.cm[2] + d * cycle.cm[0];
            mpz_mul_si(temp_0, cycle.cm[1], d);
            mpz_add(cycle.cm[3], cycle.cm[3], temp_0);
            //cycle.cm[3] = cycle.cm[3] + d * cycle.cm[1];
            mpz_mul_2exp(cycle.cm[0], cycle.cm[0], 1);
            //cycle.cm[0] = cycle.cm[0] * 2;
            mpz_mul_2exp(cycle.cm[1], cycle.cm[1], 1);
            //cycle.cm[1] = cycle.cm[1] * 2;

            mpz_mul_si(temp_0, cycle.f0_low, d);                      //按照d修改低位信息
            mpz_add(cycle.g0_low, cycle.g0_low, temp_0);
            mpz_mul_si(temp_1, cycle.f1_low, d);
            mpz_add(cycle.g1_low, cycle.g1_low, temp_1);
            mpz_mul_2exp(cycle.f0_low, cycle.f0_low, 1);
            mpz_mul_2exp(cycle.f1_low, cycle.f1_low, 1);

            mpz_mul_si(temp_0, cycle.f0_high, d);                      //按照d修改高位信息
            mpz_add(cycle.g0_high, cycle.g0_high, temp_0);
            mpz_mul_si(temp_1, cycle.f1_high, d);
            mpz_add(cycle.g1_high, cycle.g1_high, temp_1);
            mpz_mul_2exp(cycle.f0_high, cycle.f0_high, 1);
            mpz_mul_2exp(cycle.f1_high, cycle.f1_high, 1);
            
            if (cycle.powerofG > cycle.powerofF) {           //按照d修改F与G的系数信息
                for (j = 0; j < 4; j++){
                    mpz_mul_2exp(cycle.cF[j], cycle.cF[j], cycle.powerofG - cycle.powerofF);
                    //cycle.cF[j] = cycle.cF[j] * (int)pow(2, cycle.powerofG - cycle.powerofF);
                }
                cycle.powerofF = cycle.powerofG;
            } else {
                for (j = 0; j < 4; j++){
                    mpz_mul_2exp(cycle.cG[j], cycle.cG[j], cycle.powerofF - cycle.powerofG);
                    //cycle.cG[j] = cycle.cG[j] * (int)pow(2, cycle.powerofF - cycle.powerofG);
                }
                cycle.powerofG = cycle.powerofF;
            }
            for(j = 0; j < 4; j++){
                mpz_mul_si(temp_0, cycle.cF[j], d);
                mpz_add(cycle.cG[j], cycle.cG[j], temp_0);
                //cycle.cG[j] = cycle.cG[j] + d * cycle.cF[j];
            }
            cycle.powerofG++;

            mpz_mul_si(temp_0, cycle.F_low, d);                         //按照d修改F与G的低位信息
            mpz_add(cycle.G_low, cycle.G_low, temp_0);
            mpz_tdiv_q_2exp(cycle.G_low, cycle.G_low, 1);

            run = 0;    //归零游程计数器
        };
        if (dv > CHANGE ) {
            Normal(&cycle);
            dv = 0;
        }
        i++;        //一轮迭代完成，将i加一
    }

    finish = clock();

    Normal(&cycle);
    //经过所有轮计算，通过output接口输出g0与g1，最佳逼近结果即为input = g0/g1
    mpz_set(output[0], cycle.g0);
    mpz_set(output[1], cycle.g1);

    //清空占用的内存空间
    mpz_clears(temp_0, temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7, temp_8, NULL);
    clearVector(&cycle);

    //结束计时并输出总耗时
    timegap = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("i = %ld, 程序结束，总耗时%f 秒\n", i, timegap);
}

void main()
{   
    //从文件读取数据
    char input[1966001];
    FILE *p;
    //p = fopen("test-rand.txt", "r");
    p = fopen("sequence.txt", "r");
    fgets(input, 1966001, p);
    fclose(p);

    //初始化要用到的一些大数
    mpz_t output[2];
    mpz_t res, powerof2;
    mpz_inits(res, powerof2, output[0], output[1], NULL);
    mpz_ui_pow_ui(powerof2,2,1966000);  //powerof2 = 2^{1966000}

    //对input进行有理逼近，并将结果输出至output的两个分量
    Rational_Approximation_ex(input, output);

    //将答案写入另一个文件“ans.txt”
    FILE *q;
    q = fopen("ans_ex.txt", "w");
    fprintf(q, "p = ");
    mpz_out_str(q, 10, output[0]);
    fprintf(q,"\n\n");
    fprintf(q, "q = ");
    mpz_out_str(q, 10, output[1]);
    fprintf(q, "\n\n");

    //计算输出结果g0/g1 在模 2^{1966000}下的值，并以2进制写入文件。
    //如果恰是输入的倒序，则有理逼近程序计算成功，g0/g1的确是输入的一个表示
    //此时如果g0 和 g1都小于2^{1966000/2}，则g0与g1就是最佳且唯一的逼近结果
    if (mpz_sgn(output[0]) < 0) {
        mpz_add(output[0], powerof2, output[0]);
    };
    if (mpz_sgn(output[1]) < 0) {
        mpz_add(output[1], powerof2, output[1]);
    };
    mpz_invert(res, output[1], powerof2);
    mpz_mul(res, output[0], res);
    mpz_mod(res, res, powerof2);
    fprintf(q, "p/q = ");
    mpz_out_str(q, 2, res);


    //清空大数占用的内存空间，并关闭文件
    mpz_clears(res, powerof2, output[0], output[1], NULL);
    fclose(q);
}