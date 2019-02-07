#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;

in vec2 fs_Pos;
out vec4 out_Col;

const vec3 lightDirection = normalize(vec3(1.0,2.0,-1.0));
const float sceneRadius = 15.0;
const int numObjects = 10;
const float distanceThreshold = 0.001;
const float speed = 0.5;
const float roationSpeed = 1.0;

// Params definining an sdf in our scene
struct sdfParams {
    int sdfType;
    int textureType;
    vec3 center;
    float radius;
    int numCraters;
    vec3 color;
    vec3 size;
};

sdfParams sdfs[numObjects];

vec3 rotateY(vec3 point, float angle) {
    return mat3(cos(angle), 0.0, sin(angle),
                0.0,        1.0, 0.0,
                -sin(angle),0.0, cos(angle)) * point;
}

vec3 rotateX(vec3 point, float angle) {
    return mat3(1.0, 0.0,        0.0,
                0.0, cos(angle), sin(angle),
                0.0, -sin(angle),cos(angle)) * point;
}

//********************************************** Bounding Cube  Functions **************//
struct boundingCube {
   int parentIndex;
   int numObjects;
   vec3 min;
   vec3 max;
};


boundingCube bCubes1;
boundingCube[8] bCubes2;
boundingCube[64] bCubes3;
boundingCube[numObjects] bCubes4;


//This function breaks the scene up into bounding cubes
//the cubes are on three levels:
// level 1 - the whole scene
// level 2 - 8 octants of whole scene
// level 3 - 8 octants of each of the level 2 octants
void initializeBoundingCubes() {
    //initialize the big cube
    bCubes1.numObjects = 0;
    bCubes1.min = vec3(0.0, 0.0, 0.0);
    bCubes1.max = vec3(0.0, 0.0, 0.0);

    //initialize the second level cubes
    for(int i = 0; i < 8; i++) {
       bCubes2[i].numObjects = 0;
       float i2 = float(i);
       //set the minimum and max to the "lowest" corner of each cube
       float xmin1 = -sceneRadius + mod(i2, 2.0) * sceneRadius;
       float ymin1 = -sceneRadius + floor(mod(i2/2.0, 2.0)) * sceneRadius;
       float zmin1 = -sceneRadius + floor(i2/4.0) * sceneRadius;
       bCubes2[i].min = vec3(xmin1, ymin1, zmin1);
       bCubes2[i].max = vec3(xmin1, ymin1, zmin1);

       //initialize the third level cubes
       for(int j = 0; j < 8; j++) {
           int index = i*8 + j;
           float j2 = float(j);
           float xmin2 = xmin1 + mod(j2, 2.0) * sceneRadius / 2.0;
           float ymin2 = ymin1 + floor(mod(j2/2.0, 2.0)) * sceneRadius / 2.0;
           float zmin2 = zmin1 + floor(j2/4.0) * sceneRadius;
           bCubes3[index].numObjects = 0;
           bCubes3[index].min = vec3(xmin2, ymin2, zmin2);
           bCubes3[index].max = vec3(xmin2, ymin2, zmin2);
           bCubes3[index].parentIndex = i;
       }
    }
}

int getBoundingCube2Index(vec3 point) {
    //check to make sure it is in our scene
    if(abs(point.x) > sceneRadius) return - 1;
    if(abs(point.y) > sceneRadius) return - 1;
    if(abs(point.z) > sceneRadius) return - 1;

    int index = 0;
    if(point.x >= 0.0) index++;
    if(point.y >= 0.0) index += 2;
    if(point.z >= 0.0) index += 4;

    return index;
}

//get the third level bounding index where point resides
int getBoundingCube3Index(vec3 point) {
    int index = getBoundingCube2Index(point);
    if(index == -1) return -1;

    index = index * 8;
    index += 1 * int(floor(2.0 * mod(point.x, sceneRadius) / sceneRadius));
    index += 2 * int(floor(2.0 * mod(point.y, sceneRadius) / sceneRadius));
    index += 4 * int(floor(2.0 * mod(point.z, sceneRadius) / sceneRadius));

    return index;
}

boundingCube getBoundingCubeForObject(sdfParams params) {
    boundingCube bCube;
    switch(params.sdfType) {
        case 0: //sphere
        case 1: //moon
        case 2: //shere intersect
        case 3: //sphere smooth blend
            bCube.min = params.center - vec3(params.radius);
            bCube.max = params.center + vec3(params.radius);
            break;
        case 4: //sphere smooth blend
            bCube.min = params.center - params.size/2.0 - 0.1;
            bCube.max = params.center + params.size/2.0 + 0.1;
            }
    return bCube;
}

void addChildToBoundingCube(in boundingCube child, inout boundingCube parent) {
    parent.numObjects++;
    if (parent.numObjects == 1) {
        parent.min = child.min;
        parent.max = child.max;
    }
    else {
        parent.min.x = min(parent.min.x, child.min.x);
        parent.min.y = min(parent.min.y, child.min.y);
        parent.min.z = min(parent.min.z, child.min.z);

        parent.max.x = max(parent.max.x, child.max.x);
        parent.max.y = max(parent.max.y, child.max.y);
        parent.max.z = max(parent.max.z, child.max.z);
    }


}


void addObjectToBoundingCubes(int objectIndex) {
    sdfParams params = sdfs[objectIndex];
    boundingCube bc = getBoundingCubeForObject(params);

    addChildToBoundingCube(bc, bCubes1);

    int b2Index = getBoundingCube2Index(params.center);
    addChildToBoundingCube(bc, bCubes2[b2Index]);

    //set the b3 cube
    int b3Index = getBoundingCube3Index(params.center);
    addChildToBoundingCube(bc, bCubes3[b3Index]);

    //set the b4 cube
    bc.parentIndex = b3Index;
    bCubes4[objectIndex] = bc;
}


sdfParams boundingCubeToSdfParams(boundingCube bc) {
    sdfParams params;
    params.sdfType = 4;
    params.center = mix(bc.min, bc.max, 0.5);
    params.size = bc.max - bc.min;
    params.color = vec3(0.9);
    return params;
}




float random1( vec2 p , vec2 seed) {
  return fract(sin(dot(p + seed, vec2(127.1, 311.7))) * 43758.5453);
}

float random1( vec3 p , vec3 seed) {
  return fract(sin(dot(p + seed, vec3(987.654, 123.456, 531.975))) * 85734.3545);
}

vec2 random2( vec2 p , vec2 seed) {
  return fract(sin(vec2(dot(p + seed, vec2(311.7, 127.1)), dot(p + seed, vec2(269.5, 183.3)))) * 85734.3545);
}

vec2 random3( vec2 p , vec2 seed) {
  return fract(sin(vec2(dot(p + seed, vec2(311.7, 127.1)), dot(p + seed, vec2(269.5, 183.3)))) * 85734.3545);
}

vec3 random3( vec3 p, vec3 seed) {
  return fract(sin(vec3(
      dot(p + seed, vec3(311.7, 127.1, 343.2)),
      dot(p + seed, vec3(269.5, 183.3, 32.1)),
      dot(p + seed, vec3(269.5, 183.3, 432.2))
  )) * 85734.3545);
}

//############################################ SDF Manipulation Functions ################3

//function to subract one sdf defined shape from another
float sdfSubtract(float distance1, float distance2) {
    return max ( -1.0 * distance1, distance2);
}

float sdfUnion(float distance1, float distance2) {
   return min(distance1, distance2);
}


//function to find the intersection of one sdf shape with another
float sdfIntersect(float distance1, float distance2) {
    return max ( distance1, distance2 );
}

//function to find the intersection of one sdf shape with another
float sdfSmoothBlend(float a, float b) {
    float k = 0.7;
    float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
    return mix( b, a, h ) - k*h*(1.0-h);
}

//  Function to calculate the ray based upn the up, eye, ref, aspect ration and screen position
vec3 getRay(vec3 up, vec3 eye, vec3 ref, float aspect, vec2 screenPos) {
    vec3 right = normalize(cross( up - eye, up));  //right vector
    float len = length(ref - eye);   //length
    vec3 vert = up * len; //normally this would also be based upon FOV tan(FOV) but we are constraing to the box
    vec3 horiz = right * aspect * len; //normally this would also be based upon FOV tan(FOV) but we are constraining to the box
    vec3 point = ref + (screenPos.x * horiz) + screenPos.y * vert;

    //calculate the ray
    return normalize(point - eye);

}

float sphereSDF(sdfParams params, vec3 point) {
    return length(point - params.center) - params.radius;
}


vec3 sphereNormal(sdfParams params, vec3 point) {
    return normalize(point - params.center);
}



void getMoonCrater(sdfParams params, int craterIndex, out sdfParams craterParams) {

    vec3 craterPlacement = normalize(random3(vec3(craterIndex,2,3), vec3(1,2,3)) * 2.0 - 1.0);
    //rotate around
    craterPlacement = rotateY(craterPlacement, roationSpeed * u_Time / 9.0);
    float craterOffset = params.radius * (1.0 + random1(vec2(craterIndex, 2.0), vec2(2,3)) * 0.3);
    craterParams.center = params.center + (craterPlacement * craterOffset);
    craterParams.radius = params.radius * random1(vec2(craterIndex, 2.0), vec2(2,3)) * 0.4;
}


float moonSDF(sdfParams params, vec3 point) {
    sdfParams craterParams;
    float distance = sphereSDF(params, point);
    for(int i = 0; i < params.numCraters; i++) {
        getMoonCrater(params, i, craterParams);
        distance = sdfSubtract(
            sphereSDF(craterParams, point),
            distance
        );
    }
    return distance;

}

float cylindarSDF(sdfParams params, vec3 point) {
    vec3 p = point - params.center;
    vec2 h = vec2(params.radius, params.size.y);

    vec2 d = abs(vec2(length(p.xz),p.y)) - h;
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));

}

vec3 cylindarNormal(sdfParams params, vec3 point) {
     vec3 p = point - params.center;
     float h = params.size.y;
     if(point.y >  h *.5)  return vec3(0,1,0);
     if(point.y <  h *.5)  return vec3(0,-1,0);

     vec2 n1 = normalize(p.xz);
     return vec3(n1.x, 0, n1.y);

}


float discSDF(sdfParams params, vec3 point) {
     //translate points so that the disc is at origin
     sdfParams params2 = params;
     params2.size.y = params2.size.y * 1.5;
     params2.radius = 4.0;
     return sdfSubtract(
         cylindarSDF(params2, point),
         cylindarSDF(params, point)

     );
}
vec3 discNormal(sdfParams params, vec3 point) {
     vec3 p = point - params.center;
     float h = params.size.y;

     vec3 normal = cylindarNormal(params, point);

     //flip for the inside of the ring
     if(abs(p.y) < h * 0.5 && length(p.xy) < params.radius)
         return normal * -1.0;

     return normal;


}




float sphereIntersectSDF(sdfParams params, vec3 point) {
    //adjust the two sphere centers
    sdfParams params1 = params;
    params1.center.x = params.center.x + params.radius / 2.0;

    sdfParams params2 = params;
    params2.center.x = params.center.x - params.radius / 2.0;

    return sdfIntersect(
        sphereSDF(params1, point),
        sphereSDF(params2, point)
    );
}



vec3 sphereIntersectNormal(sdfParams params, vec3 point) {
    vec3 center1 = vec3(params.center.x + params.radius / 2.0, params.center.yz);
    vec3 center2 = vec3(params.center.x - params.radius / 2.0, params.center.yz);

    if(point.x > params.center.x) return normalize(point - center2);
    return normalize(point - center1);
}


float spikeSDF(vec3 center, float width, float height, vec3 direction, vec3 point) {
    if(point.y > center.y + height) return length(point - vec3(center.x, center.y + height, center.z));
    if(point.y < center.y - height) return length(point - vec3(center.x, center.y - height, center.z));

    float radius = width * (1.0 - abs(point.y - center.y)/height);
    return length(point - center);
}

float sphereSmoothBlendSDF(sdfParams params, vec3 point) {
    //create sphere above the given center
    float separation = sin(speed * u_Time/20.0) * params.radius *  1.5;
    float radius = params.radius;

    sdfParams params1 = params;
    params1.center.y += separation;
    params1.radius = radius;
    float dist1 = sphereSDF(params1, point);

    //create second spehere below the current center
    sdfParams params2 = params;
    params2.center.y -= separation;
    params2.radius = radius;
    float dist2 = sphereSDF(params2, point);

    //return dist1;

    return sdfSmoothBlend(dist1, dist2);

}

float cubeSDF(sdfParams params, vec3 point) {
     //translate points so that cube is at origin
     point -= params.center;
     vec3 d = abs(point) - params.size/2.0;
     return length(max(d, 0.0))
        + min(max(d.x,max(d.y,d.z)),0.0);
}

float rayMarch(sdfParams params, vec3 ray, int maxIterations, float maxT) {
    float t = 0.0;
    vec3 rayPos;
    float distance;
    int iterations = 0;
    while (t < maxT && iterations <= maxIterations) {

        rayPos = u_Eye + t * ray;

        //get distance from point on the ray to the object
        switch(params.sdfType) {
            case 0: distance = sphereSDF           (params, rayPos); break;
            case 1: distance = moonSDF             (params, rayPos); break;
            case 2: distance = sphereIntersectSDF  (params, rayPos); break;
            case 3: distance = sphereSmoothBlendSDF(params, rayPos); break;
            case 4: distance = cubeSDF             (params, rayPos); break;
            case 5: distance = discSDF             (params, rayPos); break;
            case 6: distance = cylindarSDF         (params, rayPos); break;
        }

        //if distance < some epsilon we are done
        if(distance < distanceThreshold) {
            return t;
        }

        t += distance;
        iterations++;
    }
    if(iterations >= maxIterations) return maxT;

    return t;
}

vec3 getNormalFromRays(sdfParams params) {
    float aspect = u_Dimensions.x / u_Dimensions.y;
    //calculate the points for 4 surrounding rays
    vec3 ray1 = getRay(u_Up, u_Eye, u_Ref, aspect, fs_Pos + vec2(-0.001,  0.0));
    vec3 ray2 = getRay(u_Up, u_Eye, u_Ref, aspect, fs_Pos + vec2( 0.001,  0.0));
    vec3 ray3 = getRay(u_Up, u_Eye, u_Ref, aspect, fs_Pos + vec2( 0.00, -0.001));
    vec3 ray4 = getRay(u_Up, u_Eye, u_Ref, aspect, fs_Pos + vec2( 0.00,  0.001));

    float t1 =  rayMarch(params, ray1, 100, 100.0);
    float t2 =  rayMarch(params, ray2, 100, 100.0);
    float t3 =  rayMarch(params, ray3, 100, 100.0);
    float t4 =  rayMarch(params, ray4, 100, 100.0);

    vec3 p1 = u_Eye + ray1 * t1;
    vec3 p2 = u_Eye + ray2 * t2;
    vec3 p3 = u_Eye + ray3 * t3;
    vec3 p4 = u_Eye + ray4 * t4;

    return normalize(cross(p4-p3, p1-p2));
}

vec3 moonNormal(sdfParams params,  vec3 point) {
   return getNormalFromRays(params);

}

vec3 cubeNormal(sdfParams params, vec3 point) {
    //return getNormalFromRays(params);
    if(abs(params.center.x + params.size.x/2.0 - point.x) <= distanceThreshold*10.0) return vec3(1.0, 0.0, 0.0);
    if(abs(params.center.x - params.size.x/2.0 - point.x) <= distanceThreshold*10.0) return vec3(-1.0, 0.0, 0.0);

    if(abs(params.center.y + params.size.y/2.0 - point.y) <= distanceThreshold) return vec3(0.0, 1.0, 0.0);
    if(abs(params.center.y - params.size.y/2.0 - point.y) <= distanceThreshold) return vec3(0.0, -1.0, 0.0);

    if(abs(params.center.z + params.size.z/2.0 - point.z) <= distanceThreshold) return vec3(0.0, 0.0, 1.0);
    if(abs(params.center.z - params.size.z/2.0 - point.z) <= distanceThreshold) return vec3(0.0, 0.0, -1.0);

    return vec3(1.0);
}


vec3 sphereSmoothBlendNormal(sdfParams params, vec3 point) {
    return getNormalFromRays(params);
}


vec3 getNormal(sdfParams params, vec3 point) {
    switch(params.sdfType) {
        case 0: return sphereNormal            (params, point);
        case 1: return moonNormal              (params, point);
        case 2: return sphereIntersectNormal   (params, point);
        case 3: return sphereSmoothBlendNormal (params, point);
        case 4: return cubeNormal              (params, point);
       // case 5: return discNormal              (params, point);
       // case 6: return cylindarNormal          (params, point);
        default: return getNormalFromRays      (params);
    }
    return vec3(0.0, 0.1, 0.0);
}


vec4 getTextureColor(sdfParams params, vec3 point) {
    vec3 normal;
    vec3 lightDirection;
    float intensity;
    switch(params.textureType) {
        ///flat lambert
        case 0:
            normal = getNormal(params, point);
            lightDirection = normalize(-1.0 * params.center);
            float intensity = dot(normal, lightDirection) * 0.9 + 0.1;
            return vec4(params.color * intensity, 1.0);
        //light from camera
        case 1:
            normal = getNormal(params, point);
            lightDirection = normalize(u_Eye - params.center);
            intensity = dot(normal, lightDirection) * 0.9 + 0.5;
            return vec4(params.color * intensity, 1.0);
             return vec4(params.color, 1.0);
             //disc
    }
    return vec4(params.color, 1.0);
}

void initSdfs() {
    //earth

    vec3 earthCenter = vec3(
        sin((u_Time+100.0)/80.0 * speed)* 18.0,
        0,
        cos((u_Time+100.0)/80.0 * speed) * 18.0
    );
    sdfs[0].center = earthCenter;
    sdfs[0].radius = 2.0;
    sdfs[0].sdfType = 0;
    sdfs[0].textureType = 0;
    sdfs[0].color = vec3(0.0, 0.0, 1.0);

    //moon
    vec3 moonCenter = earthCenter +
    vec3(
        sin((u_Time+100.0)/10.0 * speed)* 5.0,
        sin((u_Time+100.0)/10.0 * speed)* 2.0,
        cos((u_Time+100.0)/10.0 * speed) * 5.0
    );
    sdfs[1].center = moonCenter;
    sdfs[1].radius = 1.0;
    sdfs[1].sdfType = 1;
    sdfs[1].textureType = 0;
    sdfs[1].color = vec3(0.776, 0.858, 0.862);
    sdfs[1].numCraters = 20;

    //football saturn
    vec3 footbalCenter = vec3(
         sin(u_Time/120.0 * speed)* 25.0,
         0,
         cos(u_Time/120.0 * speed) * 25.0
     );
    sdfs[2].center = footbalCenter;
    sdfs[2].radius = 2.0;
    sdfs[2].sdfType = 2;
    sdfs[2].textureType = 0;
    sdfs[2].color = vec3(1.0, 0.0, 0.0);
    //football disc
    sdfs[4].center = footbalCenter;
    sdfs[4].radius = 5.0;
    sdfs[4].sdfType = 5;
    sdfs[4].textureType = 0;
    sdfs[4].color = vec3(0.3, 1.0, 1.0);
    sdfs[4].size = vec3(2.8,0.1,0.0);


    //sphere blend
    sdfs[3].center = vec3(0.0, 0.0, 0.0);
    sdfs[3].radius = 2.0;
    sdfs[3].sdfType = 3;
    sdfs[3].textureType = 1;
    sdfs[3].color = vec3(0.988, 0.992, 0.588);


    vec3 venusCenter = vec3(
        sin((u_Time+10.0)/10.0 * speed)* 8.0,
        0,
        cos((u_Time+10.0)/10.0 * speed) * 8.0
    );
    sdfs[5].center = venusCenter;
    sdfs[5].radius = 0.7;
    sdfs[5].sdfType = 0;
    sdfs[5].textureType = 0;
    sdfs[5].color = vec3(0.815, 0.670, 0.423);


}





void main() {


    float aspect = u_Dimensions.x / u_Dimensions.y;

    //calculate the ray
    vec3 ray = getRay(u_Up, u_Eye, u_Ref, aspect, fs_Pos);

    //set the default background color
    vec4 color = vec4(0.1 * (ray + vec3(1.0, 1.0, 1.0)), 1.0);


    //set up
    initSdfs();

    float maxT = 100.0;
    int maxIterations = 100;
    vec3 normal = lightDirection;
    float t;
    //get the diffuse term
    //maxT = t;

    ///here is the bounding box imnplemntation but it is slower
    //so I didn't useo

    //initializeBoundingCubes();
    for(int i=0; i<numObjects; i++) {
    // addObjectToBoundingCubes(i);
    }

    for(int i = 0; i < numObjects; i++) {
        t = rayMarch(sdfs[i], ray, maxIterations, maxT);
        if( t < maxT) {
            //get the diffuse term
            color = getTextureColor(sdfs[i], u_Eye + ray*t);
            maxT = t;
        }
    }
    //get the params to ray march highest bounding box
    sdfParams params = boundingCubeToSdfParams(bCubes1);
    t = rayMarch(params, ray, maxIterations, maxT);
    if(t < maxT) {
        for(int i = 0; i < numObjects; i++) {
            t = rayMarch(sdfs[i], ray, maxIterations, maxT);
            if( t < maxT) {
                //get the diffuse term
                color = getTextureColor(sdfs[i], u_Eye + ray*t);
                maxT = t;
            }
        }
        //loop through the next set of bounding boxes
//        for(int i = 0; i < 8; i++) {
//            if(bCubes2[i].numObjects > 0) {
//                params = boundingCubeToSdfParams(bCubes2[i]);
//                t = rayMarch(params, ray, maxIterations, maxT);
//                if(t < maxT) {
//
//                    //loop through next set of bounding boxes
//                    for(int j=0; j < 64; j++) {
//                        if(bCubes3[j].numObjects > 0 && bCubes3[j].parentIndex == i) {
//                            params = boundingCubeToSdfParams(bCubes3[j]);
//                            t = rayMarch(params, ray, maxIterations, maxT);
//                            if(t < maxT) {
//                                color = getTextureColor(params, u_Eye + ray*t);
//                            }
//                        }
//                    }
//                }
//            }
//        }
    }



    out_Col = color;
}
