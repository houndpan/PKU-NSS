#include<gmp.h>
#include<stdio.h>
#include<time.h>
#include<string.h>

void Max_mpz(mpz_t rop, mpz_t op1, mpz_t op2)
//将op1与op2中比较大的数赋值给rop
{
    if (mpz_cmp(op1, op2) != 0) {
        mpz_set(rop, op1);
    } else {
        mpz_set(rop, op2);
    };
}

int Max_cmp_pair(mpz_t op00, mpz_t op01, mpz_t op10, mpz_t op11, mpz_t max0, mpz_t max1, mpz_t temp0, mpz_t temp1) 
//比较两组数绝对值化后，最大数在哪个数组
//如果|op00|或|op01|最大，则返回1, 否则返回0
{
    mpz_abs(temp0, op00);
    mpz_abs(temp1, op01);
    Max_mpz(max0, temp0, temp1);
    mpz_abs(temp0, op10);
    mpz_abs(temp1, op11);
    Max_mpz(max1, temp0, temp1);
    if (mpz_cmp(max0, max1) != 0){
        return 1;
    } else {
        return 0;
    }
}

int Minimize(mpz_t op00, mpz_t op01, mpz_t op10, mpz_t op11, mpz_t rop, mpz_t temp_f, mpz_t temp_g, mpz_t temp00, mpz_t temp01, mpz_t temp10, mpz_t temp11, mpz_t temp1,  mpz_t temp2)
//计算使得两组数组合后绝对值最小的奇数
{  
    int d;
    if (mpz_sgn(op10) != -1) {
        if (mpz_sgn(op11) != -1){
            mpz_add(temp_f, op00, op01);
            mpz_neg(temp_f, temp_f);
            mpz_add(temp_g, op10, op11);
            mpz_fdiv_q(rop, temp_f, temp_g);
            d = mpz_get_si(rop);
        } else {
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
        d = d - 1;
    }
    mpz_set(temp00, op00);
    mpz_set(temp01, op01);
    mpz_mul_si(temp1, op10, d);
    mpz_add(temp00, temp00, temp1);
    mpz_mul_si(temp2, op11, d);
    mpz_add(temp01, temp01, temp2);
    mpz_set(temp10, temp00);
    mpz_set(temp11, temp01);
    mpz_mul_si(temp1, op10, 2);
    mpz_add(temp10, temp10, temp1);
    mpz_mul_si(temp2, op11, 2);
    mpz_add(temp11, temp11, temp2);
    if (Max_cmp_pair(temp10, temp11, temp00, temp01, temp1, temp2, temp_f, temp_g) > 0){
        return d;
    } else {
        return d + 2;
    };
}

void Rational_Approximation(char* Input, mpz_t output[2])
{
    int d;
    double v;
    int Length = strlen(Input);
    char state = Input[0];
    int i = 1;
    clock_t start, finish;
    double timegap;
    start = clock();
    while(state == '0'){        //读前若干个0
        i++;
        state = Input[i-1];
    };
    mpz_t temp_0, temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7, temp_8;           //这些是程序中可能用到的一些中间变量的存储位置，避免频繁分配空间，加速程序
    mpz_inits(temp_0, temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7, temp_8, NULL);


    mpz_t f0,f1,g0,g1,alpha;              //初始化算法的主要迭代参数
    mpz_inits(f0,f1,g0,g1,alpha, NULL);
    mpz_set_si(f0,0);
    mpz_set_si(f1,2);
    mpz_ui_pow_ui(g0,2,i-1);
    mpz_set_si(g1,1);
    mpz_ui_pow_ui(alpha, 2, i - 1);
    
    while(i < Length){
        if (i % 100000 == 0) {
            finish = clock();
            timegap = (double)(finish - start)/ CLOCKS_PER_SEC;
            v = i/100000;
            printf("i = %d, 目前耗时%f 秒\n", i, timegap / (v * v));
        }
        if (Input[i] == '1'){
            mpz_ui_pow_ui(temp_0, 2, i);
            mpz_add(alpha, alpha, temp_0);
            //gmp_printf("%Zd\n", temp_2);
        }
            mpz_mul(temp_1, alpha, g1);
            mpz_sub(temp_1,temp_1,g0);
            //gmp_printf("%Zd\n", alpha);
        if (mpz_divisible_2exp_p(temp_1, i + 1) != 0 ) {
            mpz_mul_2exp(f0,f0,1);
            mpz_mul_2exp(f1,f1,1);
        } else if (Max_cmp_pair(f0, f1, g0, g1, temp_0, temp_1, temp_2, temp_3) != 0) {
            d = Minimize(f0, f1, g0, g1, temp_0, temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7, temp_8);
            mpz_mul_si(temp_0, g0, d);
            mpz_add(f0, f0, temp_0);
            mpz_mul_si(temp_1, g1, d);
            mpz_add(f1, f1, temp_1);
            mpz_mul_2exp(g0, g0, 1);
            mpz_mul_2exp(g1, g1, 1);
            mpz_swap(g0, f0);
            mpz_swap(g1, f1);
        } else {
            d = Minimize(g0, g1, f0, f1, temp_0, temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7, temp_8);
            mpz_mul_si(temp_0, f0, d);
            mpz_add(g0, g0, temp_0);
            mpz_mul_si(temp_1, f1, d);
            mpz_add(g1, g1, temp_1);
            mpz_mul_2exp(f0, f0, 1);
            mpz_mul_2exp(f1, f1, 1);
        };
           // gmp_printf("%Zd\n", temp_2);
       // gmp_printf("%Zd,\t%Zd,\t\t%Zd,\t%Zd\n",f0, f1, g0, g1);
       //gmp_printf("i=%d,\talpha=%Zd\n",i,alpha);
		i++;
	};
    mpz_set(output[0], g0);
    mpz_set(output[1], g1);
    mpz_clears(f0,f1,g0,g1,alpha, NULL);
    mpz_clears(temp_0, temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7, temp_8, NULL);
    finish = clock();
    timegap = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("i = %d,%f \n", i, timegap);
}

int main()
{
    char input[1966000];
    FILE *p;
    p = fopen("sequence.txt", "r");
    fgets(input, 1966000, p);
    fclose(p);
    mpz_t output[2];
    mpz_t res, powerof2;
    mpz_inits(res, powerof2, output[0], output[1], NULL);
    mpz_ui_pow_ui(powerof2,2,1966000);
    Rational_Approximation(input, output);
    FILE *q;
    q = fopen("ans_old.txt", "w");
    mpz_out_str(q, 10, output[0]);
    fprintf(q,"\n\n");
    mpz_out_str(q, 10, output[1]);
    fprintf(q, "\n\n");
    if (mpz_sgn(output[0]) < 0) {
        mpz_add(output[0], powerof2, output[0]);
    };
    if (mpz_sgn(output[1]) < 0) {
        mpz_add(output[1], powerof2, output[1]);
    };mpz_invert(res, output[1], powerof2);
    mpz_mul(res, output[0], res);
    mpz_mod(res, res, powerof2);
    //mpz_mod(res, res, powerof2);
    mpz_out_str(q, 2, res);
    mpz_clears(res, powerof2, output[0], output[1], NULL);
    fclose(q);
    return 0;

}