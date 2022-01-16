class Vector {
  constructor(x, y, z) {
    this.values = [x, y, z];
  }
  vector_scale(scale) {
    let scaled_vector = new Vector(0,0,0);
    this.values.map((value, index) => {
      return scaled_vector.values[index] = value * scale;
    })
    return scaled_vector;
  }
}

function vector_addition() {
  added_vector = new Vector(0,0,0);
  for (let i = 0; i < arguments.length; i++) {
    for (let j = 0; j < 3; j++) {
      added_vector.values[j] += arguments[i].values[j];
    }
  }
  return added_vector;
}

function vector_dot_product(vector1, vector2) {
  if (vector1.length !== vector2.length) {
    return null;
  } else {
    let dot_sum = 0;
    for (let i = 0; i < vector1.values.length; i++) {
      dot_sum += vector1.values[i] * vector2.values[i];
    }
    return dot_sum;
  }
}

function vector_component_multiplication(vector1, vector2) {
  if (vector1.length !== vector2.length) {
    return null;
  } else {
    let return_vector = new Vector(0,0,0);
    for (let i = 0; i < vector1.values.length; i++) {
       return_vector.values[i] = vector1.values[i] * vector2.values[i];
    }
    return return_vector;
  }
}

// routines for creating a ray tracing scene
let scene = {};

// NEW COMMANDS FOR PART B

// create a new disk
function new_disk (x, y, z, radius, nx, ny, nz, dr, dg, db, k_ambient, k_specular, specular_pow) {
  if (!scene["disks"]) {
    scene["disks"] = [];
  }
    let disks = scene["disks"];
    disks.push({x, y, z, radius, nx, ny, nz, dr, dg, db, k_ambient, k_specular, specular_pow});
    scene["disks"] = disks;
}

// create a new area light source
function area_light (r, g, b, x, y, z, ux, uy, uz, vx, vy, vz) {
  if (!scene["arealights"]) {
    scene["arealights"] = [];
  }
  if (scene["arealights"].length <= 10) {
    let lights = scene["arealights"];
    lights.push({r, g, b, x, y, z, ux, uy, uz, vx, vy, vz});
    scene["arealights"] = lights;
  }
}

function set_sample_level (num) {
  scene["sample_level"] = num;
}

function jitter_on() {
  scene["jitter"] = true;
}

function jitter_off() {
  scene["jitter"] = false;
}


// OLD COMMANDS FROM PART A (some of which you will still need to modify)

// clear out all scene contents
function reset_scene() {
  let old_jitter = scene["jitter"];
  let old_sample = scene["sample_level"];
  if (!old_sample) {
    old_sample = 1;
  }
  if (!old_jitter) {
    old_jitter = false;
  }
  /* if (!old_sample) {
    scene = {old_jitter, old_sample, lights: [], arealights: [], sample_level: 1};
  } else {
    scene = {old_jitter, old_sample, lights: [], arealights: [], sample_level: old_sample}; 
  } */
  scene = {old_jitter, old_sample, lights: [], arealights: [], sample_level: old_sample, jitter: old_jitter};
}

// create a new point light source
function new_light (r, g, b, x, y, z) {
  if (!scene["lights"]) {
    scene["lights"] = [];
  }
  if (scene["lights"].length <= 10) {
    let lights = scene["lights"];
    lights.push({r, g, b, x, y, z});
    scene["lights"] = lights;
  }
}

// set value of ambient light source
function ambient_light (r, g, b) {
  scene["ambient"] = createVector(r, g, b);
}

// set the background color for the scene
function set_background (r, g, b) {
  scene["background"] = {r, g, b};
}

// set the field of view
function set_fov (theta) {
  scene["fov"] = theta;
}

// set the position of the virtual camera/eye
function set_eye_position (x, y, z) {
  scene["eye"] = createVector(x, y, z);
}

// set the virtual camera's viewing direction
function set_uvw(x1,y1, z1, x2, y2, z2, x3, y3, z3) {
  // should be normalized first
  let magU = Math.min(Math.sqrt(x1 * x1 + y1 * y1 + z1 * z1), 1);
  let magV = Math.min(Math.sqrt(x2 * x2 + y2 * y2 + z2 * z2), 1);
  let magW = Math.min(Math.sqrt(x3 * x3 +  y3 * y3 +  z3 * z3), 1);
  scene["u"] = createVector(x1, y1, z1).normalize();
  scene["v"] = createVector(x2, y2, z2).normalize();
  scene["w"] = createVector(x3, y3, z3).normalize();
}

// create a new sphere
function new_sphere (x, y, z, radius, dr, dg, db, k_ambient, k_specular, specular_pow) {
  if (!scene["spheres"]) {
    scene["spheres"] = [];
  }
    let spheres = scene["spheres"];
    spheres.push({x, y, z, radius, dr, dg, db, k_ambient, k_specular, specular_pow});
    scene["spheres"] = spheres;
}

// create an eye ray based on the current pixel's position
// returns an object with origin and direction (which is an object with dx, dy, dz)
function eye_ray_uvw (i, j) {
  let d = 1 / Math.tan((scene["fov"]/2) * Math.PI / 180);
  let u = -1 + 2 * i / width;
  let v = (j - height/2) * (-2 / height);
  let scaled_u = p5.Vector.mult(scene["u"], u);
  let scaled_v = p5.Vector.mult(scene["v"], v);
  let scaled_w = p5.Vector.mult(scene["w"], -d);
  let ray_direction = p5.Vector.add(p5.Vector.add(scaled_w, scaled_u), scaled_v);
  return { origin: scene["eye"], direction: ray_direction};
}

// this is the main routine for drawing your ray traced scene
function draw_scene() {

  noStroke();  // so we don't get a border when we draw a tiny rectangle

  // go through all the pixels in the image
  
  for (let y = 0; y < height; y++) {
    for (let x = 0; x < width; x++) {
      
      
      let level = scene.sample_level;
      
      /* // subdivide x and y
      let subXStart = x - .5 + (1 / (2 * scene.sample_level));
      let subXMax = x + .5 - (1 / (2 * scene.sample_level));
      let subXInc = 1 / scene.sample_level;
      
      let subYStart = y - .5 + (1 / (2 * scene.sample_level));
      let subYMax = y + .5 - (1 / (2 * scene.sample_level));
      let subYInc = 1 / scene.sample_level; */
      
      // total color sum for this pixel
      let color_at_pixel = {
        r: 0,
        g: 0,
        b: 0
      }
      
      for (let levelX = 0; levelX < scene.sample_level; levelX += 1) {
        for (let levelY = 0; levelY < scene.sample_level; levelY += 1) {
          // calculate subX for this level
          let subX = (x - .5) + (1 / (2 * scene.sample_level)) + (levelX/scene.sample_level);
          let subY = (y - .5) + (1 / (2 * scene.sample_level)) + (levelY/scene.sample_level);

          // color for this subpixel
          let color_at_subpixel = createVector(0,0,0);

          let ray = eye_ray_uvw(subX, subY);

          // find s and t
          let paramS = (subX - x) * 2;
          let paramT = (subY - y) * 2;
    
          // Figure out the pixel's color here (FOR YOU TO WRITE!!!)
          // go through scene objects and find intersections, get all t values and see which one happened first
          let shape_and_time = find_intersections(ray);
          // If intersection has been found, get color from shade function
          // If not, get color from background
          if (shape_and_time["t"]) {
            // draw color at this pixel based on sphere
            let sphere = shape_and_time["sphere"]
            let disk = shape_and_time["disk"];
            let t = shape_and_time["t"];
            if (shape_and_time["shape"] === "sphere") {
              color_at_subpixel = shadeSphere(sphere, ray, t, paramS, paramT);
            } else if (shape_and_time["shape"] === "disk") {
              color_at_subpixel = shade_disk(disk, ray, t, paramS, paramT);
            }
          } else {
            // draw background color at this pixel
            color_at_subpixel = createVector(scene["background"]["r"], scene["background"]["g"], scene["background"]["b"])
          }
          let r,g,b;  // placeholders to store the pixel's color
          r = color_at_subpixel.x;
          g = color_at_subpixel.y;
          b = color_at_subpixel.z;

          color_at_pixel.r += r;
          color_at_pixel.g += g;
          color_at_pixel.b += b;
        }
      }
      let scale_factor = 1 / (scene.sample_level * scene.sample_level);

      let r = color_at_pixel.r / (scene.sample_level * scene.sample_level);
      let g = color_at_pixel.g / (scene.sample_level * scene.sample_level);
      let b = color_at_pixel.b / (scene.sample_level * scene.sample_level);

      // set the pixel color, converting values from [0,1] into [0,255]
      fill (255 * r, 255 * g, 255 * b);
      
      rect (x, y, 1, 1);   // make a little rectangle to fill in the pixel
    }
  }
  console.log(scene);
}

function ray_intersect_sphere(sphere, ray) {
  let dx = ray.direction.x;
  let dy = ray.direction.y;
  let dz = ray.direction.z;
  let xc = sphere["x"];
  let yc = sphere["y"];
  let zc = sphere["z"];
  let xo = ray.origin.x;
  let yo = ray.origin.y;
  let zo = ray.origin.z;
  let r = sphere["radius"];
  let a = dx * dx + dy * dy + dz * dz;
  let b = ((xo - xc) * dx + (yo-yc) * dy + (zo-zc) * dz) * 2;
  let c = (xo - xc) * (xo - xc) +  (yo - yc) * (yo - yc) + (zo - zc) * (zo - zc) - r * r;
  // now solve quadratic formula
  let under_root_part = (b * b) - (4 * a * c);
  let t = 0;
  if (under_root_part > 0) {
    let t1 = (-b + Math.sqrt(under_root_part)) / (2 * a);
    let t2 = (-b - Math.sqrt(under_root_part)) / (2 * a);
    t = t1 > t2 ? t2 : t1;
  } else if (under_root_part < 0) {
    // no solutions, return negative t value
    t = -1;
  } else {
    // one solution
    t = -b / (2 * a);
  }
  return t;
}

function shadeSphere(sphere, ray, t, paramS, paramT) {
  // find point based on ray
  let initial_point = ray.origin;
  let displacement = p5.Vector.mult(ray.direction, t);
  let ray_intersect_point = p5.Vector.add(initial_point, displacement);


  // find normal
  // should be the point of intersection - the center of the sphere
  let sphere_center = createVector(sphere["x"], sphere["y"], sphere["z"]);
  let normal_vector = (p5.Vector.add(ray_intersect_point, p5.Vector.mult(sphere_center, -1))).normalize();
  


  // calculate color based on normal
  // just diffuse for now
  let diffuse = diffuse_shade(ray_intersect_point, normal_vector, createVector(sphere["dr"], sphere["dg"], sphere["db"]), sphere.k_ambient, paramS, paramT);
  if (scene.ambient) {
    let ca = new Vector(scene.ambient.x, scene.ambient.y, scene.ambient.z);
    let cr = new Vector(sphere.dr, sphere.dg, sphere.db);
    let cacr = vector_component_multiplication(ca, cr);
    let scaled = cacr.vector_scale(sphere.k_ambient);
    let ambientAddition = createVector(scaled.values[0], scaled.values[1], scaled.values[2]);
    diffuse.add(ambientAddition);
  }
  return diffuse;
}

function ray_intersect_disk(disk, ray) {

  // first find equation of plane containing disk
  let disk_normal = createVector(disk["nx"], disk["ny"], disk["nz"]).normalize();
  let disk_center = createVector(disk["x"], disk["y"], disk["z"]);
  let ax = disk_center.x * disk_normal.x;
  let by = disk_center.y * disk_normal.y;
  let cz = disk_center.z * disk_normal.z;
  let d = - (ax + by + cz);
  let plane = {a: disk_normal.x, b: disk_normal.y, c: disk_normal.z, d};

  // find point where ray intersects with plane
  let p5RayDirection = ray.direction;
  let p5RayOrigin = ray.origin;
  let x0 = p5RayOrigin.x;
  let y0 = p5RayOrigin.y;
  let z0 = p5RayOrigin.z;
  let dx = p5RayDirection.x;
  let dy = p5RayDirection.y;
  let dz = p5RayDirection.z;
  let tNumerator = - (plane.a * x0 + plane.b * y0 + plane.c * z0 + plane.d);
  let tDenominator = plane.a * dx + plane.b * dy + plane.c * dz;
  let t = tNumerator / tDenominator;
  let planeIntersect = p5.Vector.add(p5RayOrigin, p5.Vector.mult(p5RayDirection, t));

  // see if that point is inside disks radius
  // find distance between planeIntersect point and disk_center
  let distance = p5.Vector.sub(planeIntersect, disk_center).mag();
  if (distance <= disk.radius) {
    return t;
  } else {
    return null;
  }
}

function shade_disk(disk, ray, t, paramS, paramT) {
  // get normal
  let disk_normal = createVector(disk.nx, disk.ny, disk.nz).normalize();
  let displacement = p5.Vector.mult(ray.direction, t);
  
  // get intersection point
  let ray_intersect_point = p5.Vector.add(ray.origin, displacement);

  // get disk colors
  let cr = createVector(disk.dr, disk.dg, disk.db);

  // get color based on diffuse shading
  let diffuse_color = diffuse_shade(ray_intersect_point, disk_normal, cr, disk.k_ambient, paramS, paramT);
  if (scene.ambient) {
    let ca = new Vector(scene.ambient.x, scene.ambient.y, scene.ambient.z);
    let cr = new Vector(disk.dr, disk.dg, disk.db);
    let cacr = vector_component_multiplication(ca, cr);
    let scaled = cacr.vector_scale(disk.k_ambient);
    let ambientAddition = createVector(scaled.values[0], scaled.values[1], scaled.values[2]);
    diffuse_color.add(ambientAddition);
  }
  return diffuse_color;
}

function diffuse_shade(intersection_point, normal, cr, k_ambient, paramS, paramT) {
  let color_sum = createVector(0,0,0);

  // this is for point lights
  for (let i = 0; i < scene["lights"].length; i++) {
    
    // get light info
    let light = scene["lights"][i];
    let light_pos = createVector(light["x"], light["y"], light["z"]);

    //create fake intersection point shifted slightly out in normal direction
    let new_intersection_point = p5.Vector.add(intersection_point, p5.Vector.mult(normal, .00001))

    // create light vector
    // if light vector hits something BEFORE light, cl is zero this iteration
    // STILL NEED TO CHECK FOR OBJECTS AFTER LIGHT
    let light_vector_non_normal = (p5.Vector.add(light_pos, p5.Vector.mult(intersection_point, -1)));
    let light_vector = (p5.Vector.add(light_pos, p5.Vector.mult(intersection_point, -1))).normalize();
    let shadow_ray = {origin: new_intersection_point, direction: light_vector_non_normal};

    let shape_and_time_sub = find_intersections(shadow_ray);
    if (shape_and_time_sub["t"] && shape_and_time_sub["t"] <= 1) {continue} // if intersection is found, skip this light for the diffuse calculation
    
    // if light_vector intersects anything BEFORE light, make cl a zero vector
    let cl = createVector(light["r"], light["g"], light["b"]);

    // calculate N . L
    let dot = Math.max(p5.Vector.dot(light_vector, normal), 0);

    let crcl = (p5.Vector.mult(cr, cl)).mult(dot);

    color_sum.add(crcl);
  }

  for (let i = 0; i < scene["arealights"].length; i++) {
    let arealight = scene.arealights[i];
    let c = createVector(arealight.x, arealight.y, arealight.z);
    let u = createVector(arealight.ux, arealight.uy, arealight.uz);
    let v = createVector(arealight.vx, arealight.vy, arealight.vz);
    
    // do jitter if necessary
    if (scene.jitter) {
      let sJitter = ((Math.random() - .5) * 2) / (scene.sample_level * 1.2)
      let tJitter = ((Math.random() - .5) * 2) / (scene.sample_level * 1.2)

      paramS += sJitter;
      paramT += tJitter;
    }

    // find point on light based on s and t
    let su = u.mult(paramS);
    let tv = v.mult(paramT);
    let sample_point = p5.Vector.add(p5.Vector.add(c, su), tv);


    // do shadow stuff for 
    let new_intersection_point = p5.Vector.add(intersection_point, p5.Vector.mult(normal, .00001));

    let light_vector = (p5.Vector.add(sample_point, p5.Vector.mult(intersection_point, -1)));

    let shadow_ray = {origin: new_intersection_point, direction: light_vector};
    let shape_and_time_sub = find_intersections(shadow_ray);
    if (shape_and_time_sub["t"] && shape_and_time_sub["t"] <= 1) {continue} // if intersection is found, skip this light for the diffuse calculation

    // normalize after creating shadow ray
    light_vector = (p5.Vector.add(sample_point, p5.Vector.mult(intersection_point, -1))).normalize();

    let cl = createVector(arealight.r, arealight.g, arealight.b);

    let dot = Math.max(light_vector.dot(normal), 0);

    let crcl = (p5.Vector.mult(cr, cl)).mult(dot);
    console.log("cl", cl);
    console.log("cr", cr);
    console.log("crcl", crcl);
    // try something else
    let clr = arealight.r;
    let clg = arealight.g;
    let clb = arealight.b;
    
    let crr = cr.x;
    let crg = cr.y;
    let crb = cr.z;

    let compr = clr * crr * dot;
    let compg = clg * crg * dot;
    let compb = clb * crb * dot;

    let new_sum = createVector(compr, compg, compb);
    console.log(clr, clg, clb, crr, crg, crb);
    color_sum.add(new_sum);
  }

  return color_sum;
}

function find_intersections(ray) {
  let shape_and_time = {};
  if (scene.spheres) {
    for (let i = 0; i < scene["spheres"].length; i++) {
      let sphere = scene["spheres"][i];
      let t = ray_intersect_sphere(sphere, ray);

      if (t > 0 && shape_and_time["t"] > t || t > 0 && !shape_and_time["t"]) { //if old time is greater than new time(or no old time), replace it and the old sphere
        shape_and_time["t"] = t;
        shape_and_time["sphere"] = sphere;
        shape_and_time["shape"] = "sphere";
      }
    }
  }
  if (scene.disks) {
    for (let i = 0; i < scene["disks"].length; i++) {
      let disk = scene["disks"][i];
      let t = ray_intersect_disk(disk, ray);

      if (t > 0 && shape_and_time["t"] > t || t > 0 && !shape_and_time["t"]) { //if old time is greater than new time(or no old time), replace it and the old sphere
        shape_and_time["t"] = t;
        shape_and_time["disk"] = disk;
        shape_and_time["shape"] = "disk";
      }
    }
  }
  return shape_and_time;
}
