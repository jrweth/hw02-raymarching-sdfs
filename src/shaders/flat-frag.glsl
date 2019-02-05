#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;

in vec2 fs_Pos;
out vec4 out_Col;

const vec3 lightDirection = normalize(vec3(1.0,2.0,-1.0));

struct sdfParams {
    int sdfType;
    int textureType;
    vec3 center;
    float radius;
    int numCraters;
    vec3 color;
};


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

float sphereSDF(sdfParams params, vec3 point) {
    return length(point - params.center) - params.radius;
}


vec3 sphereNormal(sdfParams params, vec3 point) {
    return normalize(point - params.center);
}



void getMoonCrater(sdfParams params, int craterIndex, out sdfParams craterParams) {
    vec3 craterPlacement = normalize(vec3(-1.0, 1.0, -2.0));
    craterParams.center = params.center + (craterPlacement * params.radius);
    craterParams.radius = params.radius * 0.3;
}


float moonSDF(sdfParams params, vec3 point) {

    vec3 craterPlacement = normalize(vec3(-1.0, 1.0, -2.0));
    sdfParams craterParams;
    int craterIndex = 1;
    getMoonCrater(params, craterIndex, craterParams);

    return sdfSubtract(
        sphereSDF(craterParams, point),
        sphereSDF(params, point)
    );
}



vec3 moonNormal(sdfParams params,  vec3 point) {

    sdfParams craterParams;
    int craterIndex = 1;
    getMoonCrater(params, craterIndex, craterParams);

    if(sphereSDF(craterParams, point) < 0.0) {
        return(craterParams.center - point);
    }

    return normalize(point - params.center);
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
    sdfParams params1 = params;
    params1.center.y = params.center.y + params.radius * 0.4;
    params1.radius = params.radius * 0.4;
    float dist1 = sphereSDF(params1, point);

    //create second spehere below the current center
    sdfParams params2 = params;
    params2.center.y = params.center.y - params.radius * 0.4;
    params2.radius = params.radius * 0.4;
    float dist2 = sphereSDF(params2, point);

    //return dist1;

    return sdfSmoothBlend(dist1, dist2);

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
        }

        //if distance < some epsilon we are done
        if(distance < 0.001) {
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


vec3 sphereSmoothBlendNormal(sdfParams params, vec3 point) {
    return getNormalFromRays(params);
}




vec3 getNormal(sdfParams params, vec3 point) {
    switch(params.sdfType) {
        case 0: return sphereNormal            (params, point);
        case 1: return moonNormal              (params, point);
        case 2: return sphereIntersectNormal   (params, point);
        case 3: return sphereSmoothBlendNormal (params, point);
    }
    return vec3(0.0, 0.1, 0.0);
}


vec4 getTextureColor(sdfParams params, vec3 point) {
    switch(params.textureType) {
        ///flat lambert
        case 0:
            vec3 normal = getNormal(params, point);
            float intensity = dot(normal, lightDirection) * 0.9 + 0.1;
            return vec4(params.color * intensity, 1.0);
    }
}





void main() {


    float aspect = u_Dimensions.x / u_Dimensions.y;

    //calculate the ray
    vec3 ray = getRay(u_Up, u_Eye, u_Ref, aspect, fs_Pos);

    //set the default background color
    vec4 color = vec4(0.5 * (ray + vec3(1.0, 1.0, 1.0)), 1.0);

    //define all of our shapes
    sdfParams sdfs[4];

    //sphere
    sdfs[0].center = vec3(0.0, 0.0, 0.0);
    sdfs[0].radius = 2.0;
    sdfs[0].sdfType = 0;
    sdfs[0].textureType = 0;
    sdfs[0].color = vec3(0.6, 1.0, 0.0);

    //moon
    sdfs[1].center = vec3(4.5, 0.0, 0.0);
    sdfs[1].radius = 2.0;
    sdfs[1].sdfType = 1;
    sdfs[1].textureType = 0;
    sdfs[1].color = vec3(1.0, 0.6, 0.0);
    sdfs[1].numCraters = 1;


    //football
    sdfs[2].center = vec3(-4.5, 0.0, 0.0);
    sdfs[2].radius = 2.0;
    sdfs[2].sdfType = 2;
    sdfs[2].textureType = 0;
    sdfs[2].color = vec3(1.0, 0.0, 0.0);

    //sphere blend
    sdfs[3].center = vec3(-8.5, 0.0, 0.0);
    sdfs[3].radius = 2.0;
    sdfs[3].sdfType = 3;
    sdfs[3].textureType = 0;
    sdfs[3].color = vec3(0.0, 0.0, 1.0);

    //set up
    float maxT = 100.0;
    int maxIterations = 100;
    vec3 normal = lightDirection;
    float t;

    //ray march each sdf to determine the
    for(int i = 0; i < 4; i++) {
        t = rayMarch(sdfs[i], ray, maxIterations, maxT);

        if( t < maxT) {
            //get the diffuse term
            color = getTextureColor(sdfs[i], u_Eye + ray*t);
            maxT = t;
        }

    }


    out_Col = color;
}
