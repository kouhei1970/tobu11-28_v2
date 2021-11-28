#include "tobu.hpp"
#include "pwm_uart.hpp"
#include "ekf_sensor.hpp"
/* Private macro -------------------------------------------------------------*/
//マルチコア設定
semaphore_t sem;

float ax,ay,az,wp,wq,wr,mx,my,mz,wqa=0.0,wpa=0.0,wra=0.0,Wqa,Wpa,Wra; 
float Mn,Md;
Matrix<float, 7 ,1> xp = MatrixXf::Zero(7,1);
Matrix<float, 7 ,1> xe = MatrixXf::Zero(7,1);
Matrix<float, 7 ,1> x_sim = MatrixXf::Zero(7,1);
Matrix<float, 7 ,7> P = MatrixXf::Identity(7,7);
Matrix<float, 6 ,1> z = MatrixXf::Zero(6,1);
Matrix<float, 3, 1> omega_m = MatrixXf::Zero(3, 1);
Matrix<float, 3, 3> Q = MatrixXf::Identity(3, 3)*1;
Matrix<float, 6, 6> R = MatrixXf::Identity(6, 6)*1;
Matrix<float, 7 ,3> G;
Matrix<float, 3 ,1> beta;
volatile float kalman_time=0.0;

float dt=0.01;
short Hzcount=1;
float duty_rr,duty_rl,duty_fl,duty_fr; 
float Duty_rr,Duty_rl,Duty_fl,Duty_fr; 
float ref_e,ref_r,ref_a,err_e=0.0,err_a=0.0,err_r=0.0;
float sk_e,sk_a,sk_r;
float olderr_e,olderr_a,olderr_r;
float dk_a,dk_e,dk_r;
float se,sa,sr;
float Ti=10000;
float Td=0.0;
float h=0.01;
volatile float dmx,dmy,dmz,dxx,dxy,dxz;
float kpe=0.50,kpa=0.50,kpr=0.50;
const uint16_t Logdatanum=48000;
float logdata[Logdatanum]={0.0};
char txbuf[256];
const uint LED_PIN = 25;

void kalman(void){
   
    float ref_t,ref_phi,ref_psi,err_psi=0.0,err_phi=0.0,err_t=0.0;
    float sk_phi,sk_t,sk_psi;
    float olderr_phi,olderr_psi,olderr_t;
    float dk_phi,dk_psi,dk_t;
    float kppsi=0.01,kpphi=0.01,kpt=0.01;
    float psi,theta,phi;
    float spsi,sphi,st;
    const uint8_t DATANUM=16;
    uint32_t logcount=0;
    uint32_t printcount=0;
    float data2MID,data4MID;
    gpio_init(LED_PIN);                       //gpioを使えるようにする
    gpio_set_dir(LED_PIN, GPIO_OUT);
  //LED_PIN=0  
  while(1)
  {
    sem_acquire_blocking(&sem);
    sem_reset(&sem, 0);
    //printf("%f\n",time);
    ref_t=(Data2-data2MID)*0.523598775598299;
    ref_phi=(Data4-data4MID)*0.52398775598299;
    dt=0.01;
    omega_m <<wp,wq,wr;
    z       <<ax,ay,az,mx,my,mz;//ここに入れる
    //--Begin Extended Kalman Filter--
    ekf(xp, xe, P, z, omega_m, Q, R, G*dt, beta, dt);
    kalman_time = kalman_time + 0.01;
   
    phi=Phl(xe);
    theta=Theta(xe);
    psi=Psi(xe);

    if(Data3>=0.41){
      //phi
      olderr_phi==err_phi;
      err_phi=(ref_phi-phi);
      if (sk_phi<=3000){
        sk_phi=sk_phi+olderr_phi;
      }
      else if(-3000<=sk_phi){
        sk_phi=sk_phi+olderr_phi;
      }
      dk_phi=(err_phi-olderr_phi)*400;
      sphi=kpphi*(err_phi+1/Ti*sk_phi+Td*dk_phi);
      
      //theta
      olderr_t==err_t;
      err_t=(ref_t-theta);
      if (sk_t<=3000){
        sk_t=sk_t+olderr_t;
      }
      else if(-3000<=sk_t){
        sk_t=sk_t+olderr_t;
      }
      dk_t=(err_t-olderr_t)*400;
      st=kpt*(err_t+1/Ti*sk_t+Td*dk_t);
    }
    else{
      sk_phi=0;
      sk_t=0;
      data2MID=Data2;
      data4MID=Data4;
    }
//スティック上
  if(Data6<0.0){
    gpio_put(LED_PIN, 1);
    if(Data3>=0.41){
      logdata[logcount++]=xe(0,0);
      logdata[logcount++]=xe(1,0);
      logdata[logcount++]=xe(2,0);
      logdata[logcount++]=xe(3,0);
      logdata[logcount++]=xe(4,0);
      logdata[logcount++]=xe(5,0);
      logdata[logcount++]=xe(6,0);
      logdata[logcount++]=wp;
      logdata[logcount++]=wq;
      logdata[logcount++]=wr;
      logdata[logcount++]=ax;
      logdata[logcount++]=ay;
      logdata[logcount++]=az;
      logdata[logcount++]=mx;
      logdata[logcount++]=my;
      logdata[logcount++]=mz;

    }
  }
  //スティック下
  else if(Data6>0.0){
    gpio_put(LED_PIN, 0);
    while(logcount>=printcount){
      if(printcount<=47500){
      for (uint8_t i=0;i<DATANUM;i++){
       // printf("%10.5f",logdata[printcount+i]);

      }
      //printf("\n");
      printcount=printcount+DATANUM;
    }
   }
  }
//スティック下
  if(Data5>0.0){ 
   gpio_put(LED_PIN, 1);
   logcount=0;
   printcount=0;
}
}
}
void MAINLOOP(void)
{

  float mx1,my1,mz1;
  const uint LED_PIN = 25;
  //LED_PIN=0  
  gpio_init(25);                       //gpioを使えるようにする
  gpio_set_dir(25, GPIO_OUT);
  //e エレベータ a エルロン r ラダー t スロットル
  pwm_clear_irq(2);
  imu_mag_data_read();
  ax=   -acceleration_mg[0]*0.001*GRAV;
  ay=   -acceleration_mg[1]*0.001*GRAV;
  az=    acceleration_mg[2]*0.001*GRAV;
  wp=    angular_rate_mdps[0]*0.001*0.017453292;
  wq=    angular_rate_mdps[1]*0.001*0.017453292;
  wr=   -angular_rate_mdps[2]*0.001*0.017453292;
  dmx=  -(magnetic_field_mgauss[0]);
  dmy=   (magnetic_field_mgauss[1]);
  dmz=  -(magnetic_field_mgauss[2]);

  //回転行列
  const float rot[9]={0.65330968, 0.75327755, -0.07589064,
                     -0.75666134, 0.65302622, -0.03194321,
                      0.02549647, 0.07829232,  0.99660436};
  //中心座標
  const float center[3]={122.37559195017053, 149.0184454603531, -138.99116060635413};
  //拡大係数
  const float zoom[3]={0.003077277151877191, 0.0031893151610213463, 0.0033832794976645804};

  //回転・平行移動・拡大
  mx1 = zoom[0]*( rot[0]*dmx +rot[1]*dmy +rot[2]*dmz -center[0]);
  my1 = zoom[1]*( rot[3]*dmx +rot[4]*dmy +rot[5]*dmz -center[1]);
  mz1 = zoom[2]*( rot[6]*dmx +rot[7]*dmy +rot[8]*dmz -center[2]);
  //逆回転
  mx = rot[0]*mx1 +rot[3]*my1 +rot[6]*mz1;
  my = rot[1]*mx1 +rot[4]*my1 +rot[7]*mz1;
  mz = rot[2]*mx1 +rot[5]*my1 +rot[8]*mz1; 
  float mag_norm=sqrt(mx*mx +my*my +mz*mz);
  mx/=mag_norm;
  my/=mag_norm;
  mz/=mag_norm;


  if(Hzcount<=3){
    Hzcount=Hzcount+1;
  }
  else{
    Hzcount=1;
    sem_release(&sem);
  }

#if 1    
    //姿勢安定化
    //最大角速度 e,a 6π,r 2π
    //最大角度30°
  ref_e=Data2*18.84955592;
  ref_a=Data4*18.84955592;
  ref_r=Data1*6.283185307;

  
    //エレベータq
  olderr_e=err_e;
  err_e=(ref_e - (wq-Wqa));
  if (sk_e<=30000){
    sk_e=sk_e+olderr_e;
  }
  else if(-30000<=sk_e){
    sk_e=sk_e+olderr_e;
  }
  dk_e=(err_e-olderr_e)*400;
  se=kpe*(err_e+1/Ti*sk_e+Td*dk_e);

  //エルロンp

  olderr_a=err_a;
  err_a=(ref_a - (wp-Wpa));
  if (sk_a<=30000){
    sk_a=sk_a+olderr_a;
  }
  else if(-30000<=sk_a){
    sk_a=sk_a+olderr_a;
  }
  dk_a=(err_a-olderr_a)*400;
  sa=kpa*(err_a+1/Ti*sk_a+Td*dk_a);

  //ラダーr

  olderr_r=err_r;
  err_r=(ref_r - (wr-Wra));
  if (sk_r<=30000){
    sk_r=sk_r+olderr_r;
  }
  else if(-30000<=sk_r){
    sk_r=sk_r+olderr_r;
  }
  dk_r=(err_r-olderr_r)*400;
  sr=kpr*(err_r+1/Ti*sk_r+Td*dk_r);


  Duty_fr=Data3+(se-sa+sr)*0.25;
  Duty_fl=Data3+(se+sa-sr)*0.25;
  Duty_rr=Data3+(-se-sa-sr)*0.25;         
  Duty_rl=Data3+(-se+sa+sr)*0.25;
#endif
#if 0   
  Duty_fr=Data3;
  Duty_fl=Data3;
  Duty_rr=Data3;         
  Duty_rl=Data3;
#endif
  tight_loop_contents();




  duty_rr=(float)(DUTYMAX-DUTYMIN)*Duty_rr+DUTYMIN;
  duty_fr=(float)(DUTYMAX-DUTYMIN)*Duty_fr+DUTYMIN+96;
  duty_rl=(float)(DUTYMAX-DUTYMIN)*Duty_rl+DUTYMIN+119;
  duty_fl=(float)(DUTYMAX-DUTYMIN)*Duty_fl+DUTYMIN;
  //if (duty_rr<DUTYMIN+100)duty_rr=DUTYMIN+5.0;
  //if (duty_rl<DUTYMIN+100)duty_rl=DUTYMIN+5.0;
  //if (duty_fr<DUTYMIN+100)duty_fr=DUTYMIN+5.0;
  //if (duty_fl<DUTYMIN+100)duty_fl=DUTYMIN+5.0;




  if (duty_rr>DUTYMAX-50.0)duty_rr=DUTYMAX-50.0;
  if (duty_rr<DUTYMIN+5.0)duty_rr=DUTYMIN+5.0;
  if (duty_fr>DUTYMAX-50.0)duty_fr=DUTYMAX-50.0;
  if (duty_fr<DUTYMIN+5.0)duty_fr=DUTYMIN+5.0;
  if (duty_rl>DUTYMAX-50.0)duty_rl=DUTYMAX-50.0;
  if (duty_rl<DUTYMIN+5.0)duty_rl=DUTYMIN+5.0;
  if (duty_fl>DUTYMAX-50.0)duty_fl=DUTYMAX-50.0;
  if (duty_fl<DUTYMIN+5.0)duty_fl=DUTYMIN+5.0;

  if(Data3<0.05)
  {
    pwm_set_chan_level(slice_num[0], PWM_CHAN_A, DUTYMIN);
    pwm_set_chan_level(slice_num[0], PWM_CHAN_B, DUTYMIN);
    pwm_set_chan_level(slice_num[1], PWM_CHAN_A, DUTYMIN);
    pwm_set_chan_level(slice_num[1], PWM_CHAN_B, DUTYMIN);
    sk_r=0.0;
    sk_a=0.0;
    sk_e=0.0;
    err_r=0.0;
    err_a=0.0;
    err_e=0.0;


  }
  else
  {
    pwm_set_chan_level(slice_num[0], PWM_CHAN_A, duty_rr);
    pwm_set_chan_level(slice_num[0], PWM_CHAN_B, duty_fr);
    pwm_set_chan_level(slice_num[1], PWM_CHAN_A, duty_rl);
    pwm_set_chan_level(slice_num[1], PWM_CHAN_B, duty_fl);
  }
  //printf("%04f %04f %04f %04f %04f %04f \n",Olddata[0],Olddata[1],Olddata[2],Olddata[3],Olddata[4],Olddata[5]);    
  //e_time=time_us_64();
  //d_time=e_time-s_time;
  //printf("右前%04f   左前%04f   右後ろ%04f   左後ろ%04f\n",duty_fr,duty_fl,duty_rr,duty_rl);
  //printf("%04f \n",t);
  //t=t+1;
  
}
int main(void)
{
  float w,f=1;
  const uint LED_PIN = 25;          //LED_PIN=0  
  gpio_init(25);                       //gpioを使えるようにする
  gpio_set_dir(25, GPIO_OUT);
  xe << 1.0, 0.0, 0.0, 0.0, -0.078, 0.0016, 0.00063;
  xp =xe;

  G <<  0.0,0.0,0.0, 
        0.0,0.0,0.0, 
        0.0,0.0,0.0, 
        0.0,0.0,0.0, 
        1.0,0.0,0.0, 
        0.0,1.0,0.0, 
        0.0,0.0,1.0;

  beta << 0.00, 0.00, 0.00;

  P <<  1,0,0,0,0,0,0,  
        0,1,0,0,0,0,0,
        0,0,1,0,0,0,0,  
        0,0,0,1,0,0,0, 
        0,0,0,0,1,0,0,  
        0,0,0,0,0,1,0,  
        0,0,0,0,0,0,1;

  Q << 1.0e-6,0,0,
       0,1.0e-6,0,
       0,0,1.0e-6;

  R << 10.0e-5,0,0,0,0,0,
       0,7.5e-5,0,0,0,0,
       0,0,11.7e-5,0,0,0,
       0,0,0,36.4e-6,0,0,
       0,0,0,0,42.5e-6,0,
       0,0,0,0,0,90.5e-6;


  //gpio_put(LED_PIN, 1);
  stdio_init_all();
  //sleep_ms(1000);
  imu_mag_init();
 // gpio_put(LED_PIN, 0);
  serial_settei();
 // gpio_put(LED_PIN, 1);
  

  /* for (i=0;i<waittime;i++)
  {
    printf("#Please wait %d[s] ! \n",waittime-i);
    sleep_ms(1000);
  }
  printf("#Start Kalman Filter\n");
*/
  wpa=0.0;
  wqa=0.0;
  wra=0.0;

  
  while(f<=400){
    float mx1,my1,mz1;
    imu_mag_data_read();
    wp=    angular_rate_mdps[0]*0.001*0.017453292;
    wq=    angular_rate_mdps[1]*0.001*0.017453292;
    wr=   -angular_rate_mdps[2]*0.001*0.017453292;
    dmx=  -(magnetic_field_mgauss[0]);
    dmy=   (magnetic_field_mgauss[1]);
    dmz=  -(magnetic_field_mgauss[2]);

    //回転行列
    const float rot[9]={0.65330968, 0.75327755, -0.07589064,
                       -0.75666134, 0.65302622, -0.03194321,
                        0.02549647, 0.07829232,  0.99660436};
    //中心座標
    const float center[3]={122.37559195017053, 149.0184454603531, -138.99116060635413};
    //拡大係数
    const float zoom[3]={0.003077277151877191, 0.0031893151610213463, 0.0033832794976645804};

    //回転・平行移動・拡大
    mx1 = zoom[0]*( rot[0]*dmx +rot[1]*dmy +rot[2]*dmz -center[0]);
    my1 = zoom[1]*( rot[3]*dmx +rot[4]*dmy +rot[5]*dmz -center[1]);
    mz1 = zoom[2]*( rot[6]*dmx +rot[7]*dmy +rot[8]*dmz -center[2]);
    //逆回転
    mx = rot[0]*mx1 +rot[3]*my1 +rot[6]*mz1;
    my = rot[1]*mx1 +rot[4]*my1 +rot[7]*mz1;
    mz = rot[2]*mx1 +rot[5]*my1 +rot[8]*mz1; 
    float mag_norm=sqrt(mx*mx +my*my +mz*mz);
    mx/=mag_norm;
    my/=mag_norm;
    mz/=mag_norm;
    
    wqa=wq+wqa;
    wpa=wp+wpa;
    wra=wr+wra;
    Mn=Mn+mx;
    Md=Md+mz;
    f=f+1;
    sleep_us(2500);
  }
  Wqa=wqa/400;
  Wpa=wpa/400;
  Wra=wra/400;
  Mn=Mn/400;
  Md=Md/400;

  printf("%f,%f\n",Mn,Md);

 // gpio_put(LED_PIN, 1);
  sem_init(&sem, 0, 1);
  multicore_launch_core1(kalman);
  
  
  pwm_settei();
  
  float old_time=-1.0;
  char txbuf[256];

  while(1)
  {
   /* if(kalman_time>old_time)
    {   
      sprintf(txbuf,"%8.2f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f "
             "%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",
             kalman_time,xe(0,0),xe(1,0),xe(2,0),xe(3,0),xe(4,0),xe(5,0),xe(6,0),
             wp,wq,wr,ax,ay,az,mx,my,mz);
      printf("%s",txbuf);
      //printf("%8.2f %9.5f %9.5f %9.5f\n",kalman_time,dmx,dmy,dmz);
      old_time=kalman_time;
    }*/
  }
    //gpio_put(LED_PIN, 1);
 
}

