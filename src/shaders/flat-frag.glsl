#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;

in vec2 fs_Pos;
out vec4 out_Col;


float circleSDF(vec3 center, float radius, vec3 point) {
    return length(point - center) - radius;
}

float rayMarchCircle(vec3 center, float radius, vec3 ray, int maxIterations, float maxT) {
    float t = 0.0;
    float distance;
    int iterations = 0;
    while (t < maxT && iterations <= maxIterations) {
        //get distance from point on the ray to the object
        distance = circleSDF(center, radius, u_Eye + t * ray);

        //if distance < some epsilon we are done
        if(circleSDF(center, radius, u_Eye + t * ray) < 0.01) {
            return t;
        }

        t += distance;
        iterations++;
    }
    if(iterations >= maxIterations) return maxT;

    return t;
}


void main() {


  float aspect = u_Dimensions.x / u_Dimensions.y;

  vec3 U = u_Up;
  vec3 R = normalize(cross( u_Ref - u_Eye, u_Up));
  float len = length(u_Ref - u_Eye);
  vec3 V = u_Up * len; //normally this would also be based upon FOV tan(FOV) but we are constraing to the box
  vec3 H = R * aspect * len; //normally this would also be based upon FOV tan(FOV) but we are constraining to the box

  vec3 p = u_Ref + (fs_Pos.x * H) + fs_Pos.y * V;

  vec3 ray = normalize(p - u_Eye);


  //ray march the circle
  vec3 center = vec3(-5,0,1);
  float radius = 2.0;
  float maxT = 100.0;
  int maxIterations = 100;
  float t = rayMarchCircle(center, radius, ray, maxIterations, maxT);

  if( t < maxT) {
      out_Col = vec4(0,0,0,1.0);
  }
  else {
      out_Col = vec4(0.5 * (ray + vec3(1.0, 1.0, 1.0)), 1.0);
  }

  float dist = circleSDF(center, radius, u_Eye);
  //out_Col = vec4(0.0, 0.0, t, 1.0);

}
