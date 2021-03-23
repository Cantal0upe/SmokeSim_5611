
final int WIDTH = 200;
final int HEIGHT = WIDTH;


//interesting visc and diff values:
//diff 0, visc 0.003 (or higher)
//diff 0.00001, visc 0
float diff = 0.000001;
float visc = 0.0;

int mouseX_prev = 0;
int mouseY_prev = 0;

int mouseX_curr = 0;
int mouseY_curr = 0;

float[][] rho = new float[HEIGHT+2][WIDTH+2];
float[][] u   = new float[HEIGHT+2][WIDTH+2];
float[][] v   = new float[HEIGHT+2][WIDTH+2];

float[][] rho_prev = new float[HEIGHT+2][WIDTH+2];
float[][] u_prev   = new float[HEIGHT+2][WIDTH+2];
float[][] v_prev   = new float[HEIGHT+2][WIDTH+2];

int effect = 5;

int count = 0;
boolean isPressed = false;
boolean isLive = true;

int cursorX = WIDTH/2 +1;
int cursorY = WIDTH/2 + 1;
int prevX = WIDTH/2 + 1;
int prevY = WIDTH/2 + 1;
float h,s,b;
void settings(){
  size(WIDTH*2, HEIGHT*2);
}

void setup(){
  colorMode(HSB);
}



void draw(){
    background(0);
    
    float dt = 1 / (frameRate);
    if(isLive == true){
        
            
            
            float uV = random(10)*(random(2)-1);
            float vV = random(10)*(random(2)-1);
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    rho[cursorY+i][cursorX+j] += 0.2;
                    u[cursorY+i][cursorX+j] = uV;
                    v[cursorY+i][cursorX+j] = vV;
                }
            }
        
    }
    mouseX_prev = mouseX_curr;
    mouseY_prev = mouseY_curr;
    
    mouseX_curr = mouseX;
    mouseY_curr = mouseY;
    
    //if(mouseX != 0) v[mouseY/2][mouseX/2] -= 3*(mouseX_prev - mouseX);
    //if(mouseY != 0) u[mouseY/2][mouseX/2] -= 3*(mouseY_prev - mouseY);
    
    float velX = (mouseX - mouseX_prev);
    float velY = (mouseY - mouseY_prev);
    
    
    if(mouseX != 0 && mouseX_prev != 0){
      int[] lin = generateLine(mouseX_prev/2, mouseY_prev/2, mouseX/2, mouseY/2);
      for(int i=0; i<lin.length/2; i++){
        int x = max(lin[2*i], 0);
        x = min(x, WIDTH+1);
        int y = max(lin[2*i+1], 0);
        y = min(y, HEIGHT+1);
        v[y][x] += velX;
        u[y][x] += velY;
      }
    }
    vel_step(u, v, u_prev, v_prev, visc, dt);
    dens_step(rho, rho_prev, u, v, diff, dt);
    
    
    
    
    
    pushMatrix();
    noStroke();
    for(int i=0; i<HEIGHT; i++){
        for(int j=1; j<WIDTH; j++){
            pushMatrix();
            noStroke();
            switch(effect){
              case 1: //"NOXIOUS GAS"
                h = (((rho[i+1][j+1]*255)+50))%360;
                s = 200;
                b = (-cos((rho[i+1][j+1]*255)*PI/30)+1)*255/2;
                diff = 0.000001;
                visc = 0.0;
                break;
              case 2: //"OIL SLICK"
                h = (((rho[i+1][j+1]*255)+50)*255)%360;
                s = 200;
                b = (-cos((rho[i+1][j+1]*255)*PI/30)+1)*255/2;
                diff = 0.002; 
                visc = 0.0;
                break;
              case 3: //"PSYCHEDELIC NIGHTMARE"
                h = (((rho[i+1][j+1]*255)+50)*255)%360;
                s = 200;
                b = (-cos((rho[i+1][j+1]*255)*PI/30)+1)*255/2;
                diff = 0.000001;
                visc = 0.0;
                break;
              case 4: //"SMOKE"
                h = (((rho[i+1][j+1]*255)+50))%360;
                s = 0;
                b = (-cos((rho[i+1][j+1]*255)*PI/60)+1)*255/2;
                diff = 0.000001;
                visc = 0.0;
                break;
              case 5: //FIRE
                float d = rho[i+1][j+1] * 255*3;
                h = abs((d/4)%60);
                //h = 0;
                //s = 255-2*((d+120)%360);
                s = (sin((d-30)*PI/180)+0.5)*255;
                b = (sin((d-30)*PI/180)+0.5)*255;
                diff = 0.00003;
                visc = 0.000;
                break;
              default:
                h = (((rho[i+1][j+1]*255)+50))%360;
                s = 0;
                b = (-cos((rho[i+1][j+1]*255)*PI/60)+1)*255/2;
                diff = 0.000001;
                visc = 0.0;
            }
            
            fill(h,s,b);
            rect(j*2-1, i*2, 2, 2);
            popMatrix();
        }
    }
    popMatrix();
    surface.setTitle("FPS: " + (int)frameRate);
}


void mouseDragged(){
    count +=1;
    cursorX = mouseX/2;
    cursorY = mouseY/2;
    if(count == 3){
        prevX = pmouseX/2;
        prevY = pmouseY/2;
        count = 0;
    }
}
void mousePressed(){
    isPressed = true;
}
void mouseReleased(){
    isPressed = false;
}
void keyPressed(){
  switch(key){
    case ' ':
        if(isLive == true){
            isLive = false;
        }
        else{
            isLive = true;
        }
        break;
      
    
    case '1':
       effect = 1;
       for(int i = 0; i <WIDTH+2; i++){
          for(int j = 0; j < HEIGHT+2; j++){
            rho[i][j] = 0;
            rho_prev[i][j] = 0;
          }
       }
       break;
    case '2':
       effect = 2; 
       for(int i = 0; i <WIDTH+2; i++){
          for(int j = 0; j < HEIGHT+2; j++){
            rho[i][j] = 0;
            rho_prev[i][j] = 0;
          }
       }
       break;
    case '3':
       effect = 3; 
       for(int i = 0; i <WIDTH+2; i++){
          for(int j = 0; j < HEIGHT+2; j++){
            rho[i][j] = 0;
            rho_prev[i][j] = 0;
          }
       }
       break;
     case '4':
       effect = 4; 
       for(int i = 0; i <WIDTH+2; i++){
          for(int j = 0; j < HEIGHT+2; j++){
            rho[i][j] = 0;
            rho_prev[i][j] = 0;
          }
       }
       break;
     case '5':
       effect = 5; 
       for(int i = 0; i <WIDTH+2; i++){
          for(int j = 0; j < HEIGHT+2; j++){
            rho[i][j] = 0;
            rho_prev[i][j] = 0;
          }
       }
  }
}

int[] generateLine(int x1, int y1, int x2, int y2){
  int [] out;
  int temp;
  if(abs(y2 - y1) > abs(x2 - x1)){
    out = generateLine(y1, x1, y2, x2);
    for(int i=0; i<out.length/2; i++){
      temp = out[2*i];
      out[2*i] = out[2*i+1];
      out[2*i+1] = temp;
    }
    return out;
  }
  
  if(x2 < x1){
    temp = x1;
    x1 = x2;
    x2 = temp;
    temp = y1;
    y1 = y2;
    y2 = temp;
  }
  int dY = 1;
  if(y2 < y1){
    dY = -1;
  }
  out = new int[2*(x2 - x1)];
  int m = dY * (2 * (y2 - y1));
  int slope_error = m -  (x2 - x1)*dY;
  int i=0;
  if(out.length == 0)return out;
  for(int x = x1, y = y1; x < x2; x++){
    out[i] = x;
    out[i+1] = y;
    slope_error += m;
    
    if(slope_error >= 0){
      y += dY;
      slope_error -= 2 * (x2 - x1);
    }
    
    i+=2;
  }
  return out;
}

void vel_step(float[][] u, float[][] v, float[][] u0, float[][] v0, float visc, float dt){
    addSource(u, u0, dt);
    addSource(v, v0, dt);
    diffuse(1, u0, u, visc, dt); 
    diffuse(2, v0, v, visc, dt);
    project(u0, v0, u, v);
    advect(1, u, u0, u0, v0, dt);
    advect(2, v, v0, u0, v0, dt);
    project(u, v, u0, v0);
}

void dens_step(float[][] x, float[][] x0, float[][] u, float[][] v, float diff, float dt){
    diffuse(0, x0, x, diff, dt);
    advect(0, x, x0, u, v, dt);
}

void copy(float[][] a, float[][] b){
  for(int i = 0; i<b.length; i++){
    for(int j = 0; j<b[0].length; j++){
      a[i][j] = b[i][j];
    }
  }
}

//x and s must be the same size
void addSource(float[][] x, float[][] s, float dt){
    for(int i=0; i<HEIGHT+2; i++){
        for(int j=0; j<WIDTH+2; j++){
            x[i][j] += s[i][j]*dt;
        }
    }
}


void diffuse(int b, float[][] x, float[][] x0, float diff, float dt ){
    float a  = dt* diff * HEIGHT * WIDTH;
    
    for(int k = 0; k < 20; k++){
        for(int i = 1; i <= HEIGHT; i++){
            for(int j = 1; j <= WIDTH; j++){
                x[i][j] = (x0[i][j] + a*(x[i-1][j]+x[i+1][j]+x[i][j-1]+x[i][j+1]))/(1+4*a);
            }
        }
        set_bnd(b, x);
    }
}

void advect (int b, float[][] d, float[][] d0, float[][] u, float[][] v, float dt )
{
  int N = WIDTH;
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;
    dt0 = dt*N;
    for ( i=1 ; i<=N ; i++ ) {
        for ( j=1 ; j<=N ; j++ ) {
            x = i-dt0*u[i][j];
            y = j-dt0*v[i][j];
            
            if (x<0.5) x=0.5;
            if (x>N+0.5) x=N+ 0.5;
            i0=(int)x;
            i1=i0+1;
            
            if (y<0.5) y=0.5;
            if (y>N+0.5) y=N+ 0.5;
            j0=(int)y; 
            j1=j0+1;
            
            s1 = x-i0;
            s0 = 1.0-s1;
            
            t1 = y-j0;
            t0 = 1-t1;
            
            
            d[i][j] = 
              s0*(t0*d0[i0][j0] + t1 * d0[i0][j1])+
              s1*(t0*d0[i1][j0] + t1 * d0[i1][j1]);
        }
    }
    set_bnd (b, d );
}





void swap(float[][] a, float[][] b){
    float [][] tmp = a;
    a = b;
    b = tmp;
}


void project ( float[][] u, float[][] v, float[][] p, float[][] div )
{
    int N = WIDTH;
    int i, j, k;
    float h;
    h = 1.0/N;
    for ( i=1 ; i<=N ; i++ ) {
        for ( j=1 ; j<=N ; j++ ) {
            div[i][j] = -0.5*h*(u[i+1][j]-u[i-1][j]+
            v[i][j+1]-v[i][j-1]);
            p[i][j] = 0;
        }
    }
    set_bnd (0, div ); set_bnd (0, p );
    for ( k=0 ; k<20 ; k++ ) {
        for ( i=1 ; i<=N ; i++ ) {
            for ( j=1 ; j<=N ; j++ ) {
                p[i][j] = (div[i][j]+p[i-1][j]+p[i+1][j]+
                p[i][j-1]+p[i][j+1])/4;
            }
        }
        set_bnd ( 0, p );
    }
    for ( i=1 ; i<=N ; i++ ) {
        for ( j=1 ; j<=N ; j++ ) {
            u[i][j] -= 0.5*(p[i+1][j]-p[i-1][j])/h;
            v[i][j] -= 0.5*(p[i][j+1]-p[i][j-1])/h;
        }
    }
    set_bnd ( 1, u ); set_bnd ( 2, v );
}



void set_bnd(int b, float[][] x){
    for(int i=0; i<WIDTH+2; i++)
        x[0       ][i] = ((b==1)?-1:1) * x[1][i];
    
    for(int i=0; i<WIDTH+2; i++)
        x[HEIGHT+1][i] = ((b==1)?-1:1) * x[HEIGHT][i];
    
    for(int i=0; i<HEIGHT+2; i++)
        x[i][0      ] = (b==2 ? -1 : 1) * x[i][1];
    
    for(int i=0; i<HEIGHT+2; i++)
        x[i][WIDTH+1] = (b==2 ? -1 : 1) * x[i][WIDTH];
    
    
    x[0       ][0      ] = 0.5 * (x[1     ][0      ] + x[0       ][1    ]);
    x[0       ][WIDTH+1] = 0.5 * (x[1     ][WIDTH+1] + x[0       ][WIDTH]);
    x[HEIGHT+1][0      ] = 0.5 * (x[HEIGHT][0      ] + x[HEIGHT+1][1    ]);
    x[HEIGHT+1][WIDTH+1] = 0.5 * (x[HEIGHT][WIDTH+1] + x[HEIGHT+1][WIDTH]);
    
}
