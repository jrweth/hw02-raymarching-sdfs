#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;

in vec2 fs_Pos;
out vec4 out_Col;

void main() {


  float aspect = u_Dimensions.x / u_Dimensions.y;

  vec3 U = u_Up;
  vec3 R = normalize(cross( u_Ref - u_Eye, u_Up));
  float len = length(u_Ref - u_Eye);
  vec3 V = u_Up * len; //normally this would also be based upon FOV tan(FOV) but we are constraing to the box
  vec3 H = R * aspect * len; //normally this would also be based upon FOV tan(FOV) but we are constraining to the box

  vec3 p = u_Ref + (fs_Pos.x * H) + fs_Pos.y * V;

  vec3 ray = normalize(p - u_Eye);

  out_Col = vec4(0.5 * (ray + vec3(1.0, 1.0, 1.0)), 1.0);
}
