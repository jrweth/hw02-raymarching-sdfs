#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;

in vec2 fs_Pos;
out vec4 out_Col;

const vec3 lightDirection = normalize(vec3(1.0,1.0,-1.0));


//float random1( vec2 p , vec2 seed) {
//  return fract(sin(dot(p + seed, vec2(127.1, 311.7))) * 43758.5453);
//}
//
//float random1( vec3 p , vec3 seed) {
//  return fract(sin(dot(p + seed, vec3(987.654, 123.456, 531.975))) * 85734.3545);
//}
//
//vec2 random2( vec2 p , vec2 seed) {
//  return fract(sin(vec2(dot(p + seed, vec2(311.7, 127.1)), dot(p + seed, vec2(269.5, 183.3)))) * 85734.3545);
//}
//
//vec2 random3( vec2 p , vec2 seed) {
//  return fract(sin(vec2(dot(p + seed, vec2(311.7, 127.1)), dot(p + seed, vec2(269.5, 183.3)))) * 85734.3545);
//}
//
//vec2 random3( vec3 p, vec3 seed) {
//  return fract(sin(vec3(
//      dot(p + seed, vec3(311.7, 127.1, 343.2)),
//      dot(p + seed, vec3(269.5, 183.3, 32.1)),
//      dot(p + seed, vec3(269.5, 183.3, 432.2))
//  )) * 85734.3545);
//}


float sphereSDF(vec3 center, float radius, vec3 point) {
    return length(point - center) - radius;
}



float rayMarchSphere(vec3 center, float radius, vec3 ray, int maxIterations, float maxT) {
    float t = 0.0;
    float distance;
    int iterations = 0;
    while (t < maxT && iterations <= maxIterations) {
        //get distance from point on the ray to the object
        distance = sphereSDF(center, radius, u_Eye + t * ray);

        //if distance < some epsilon we are done
        if(distance < 0.01) {
            return t;
        }

        t += distance;
        iterations++;
    }
    if(iterations >= maxIterations) return maxT;

    return t;
}

vec3 getSphereNormal(vec3 center, vec3 point) {
    return normalize(point - center);
}


void getMoonCrater(vec3 moonCenter, float moonRadius, int craterIndex, out vec3 craterCenter, out float craterRadius) {
    vec3 craterPlacement = normalize(vec3(-1.0, 1.0, -2.0));
    craterCenter = moonCenter + (craterPlacement * moonRadius);
    craterRadius = moonRadius * 0.3;
}


float moonSDF(vec3 center, float radius, int numCraters, vec3 point) {

    vec3 craterPlacement = normalize(vec3(-1.0, 1.0, -2.0));
    vec3 craterCenter;// = center + (craterPlacement * radius);
    float craterRadius;// = radius * 0.5;
    int craterIndex = 1;
    getMoonCrater(center, radius, craterIndex, craterCenter, craterRadius);

    return max (
        -1.0 * sphereSDF(craterCenter, craterRadius, point),
        sphereSDF(center, radius, point)
    );

}

float rayMarchMoon(
    vec3 center,
    float radius,
    int numCraters,
    vec3 ray,
    int maxIterations,
    float maxT
) {
    float t = 0.0;
    float distance;
    int iterations = 0;
    while (t < maxT && iterations <= maxIterations) {
        //get distance from point on the ray to the object
        distance = moonSDF(center, radius, numCraters, u_Eye + t * ray);

        //if distance < some epsilon we are done
        if(distance < 0.01) {
            return t;
        }

        t += distance;
        iterations++;
    }
    if(iterations >= maxIterations) return maxT;

    return t;
}


vec3 getMoonNormal(vec3 center, float radius, int numCraters, vec3 point) {

    vec3 craterCenter;
    float craterRadius;
    int craterIndex = 1;
    getMoonCrater(center, radius, craterIndex, craterCenter, craterRadius);

    if(sphereSDF(craterCenter, craterRadius, point) < 0.0) {
        return(craterCenter - point);
    }

    return normalize(point - center);
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


  //ray march the sphere
  vec3 center = vec3(0,0,0);
  float radius = 2.0;
  int numCraters = 1;
  float maxT = 100.0;
  int maxIterations = 100;
  vec3 color = vec3(0.5, 0.5, 1.0);
  float t = rayMarchMoon(center, radius, numCraters, ray, maxIterations, maxT);

  if( t < maxT) {
      //get the diffuse term
      vec3 normal = getMoonNormal(center, radius, numCraters, u_Eye + ray*t);
      //get the lambert intesity based upon the normal
      float intensity = dot(normal, lightDirection) * 0.9 + 0.1;
      out_Col = vec4(color * intensity, 1.0);
  }
  else {
      out_Col = vec4(0.5 * (ray + vec3(1.0, 1.0, 1.0)), 1.0);
  }




  //float dist = sphereSDF(center, radius, u_Eye);
  //out_Col = vec4(0.0, 0.0, t, 1.0);

}
