let eps = 0.0001;
let timestep = 10;
let last_time = null;

/*
================
LINEAR ALGEBRA
================
*/

function matinv(M) {
    console.log(M.toString());
    let n = M.length;
    let Maug = []
    for (let row = 0; row < n; row++) {
        Maug.push([...M[row]]);
        for (let aug_col = 0; aug_col < n; aug_col++) {
            Maug[Maug.length-1].push((aug_col == row)? 1 : 0);
        }
    }
    for (let col = 0; col < n; col++) {
        let best = 0;
        let bestRow = undefined;
        for (let row = col; row < n; row++) {
            if (Math.abs(Maug[row][col]) > Math.abs(best)) {
                bestRow = row;
                best = Maug[row][col];
            }
        }
        for (let row = 0; row < n; row++) {
            if (bestRow != row) {
                let factor = Maug[row][col] / best;
                for (let subCol = 0; subCol < 2*n; subCol++) {
                    Maug[row][subCol] -= factor * Maug[bestRow][subCol];
                }
            }
        }
        for (let subCol = 0; subCol < 2*n; subCol++) {
            Maug[bestRow][subCol] /= best;
        }
        if (bestRow != col) {
            let temp = Maug[col];
            Maug[col] = Maug[bestRow];
            Maug[bestRow] = temp;
        }
    }
    let Minv = [];
    for (let row = 0; row < n; row++) {
        Minv.push([...Maug[row]]);
        Minv[row].splice(0, n);
    }
    return Minv;
}
function mattrans(M) {
    let MT = [];
    for (let i = 0; i < M[0].length; i++) {
        let row = [];
        for (let j = 0; j < M.length; j++) {
            row.push(M[j][i]);
        }
        MT.push(row);
    }
    return MT;
}
function vecdot(v1, v2) {
    let out = 0;
    for (let i = 0; i < v1.length; i++) {
        out += v1[i] * v2[i];
    }
    return out;
}
function matvecmul(M, v) {
    let Mv = [];
    for (let row of M) {
        Mv.push(vecdot(row, v));
    }
    return Mv;
}

/*
================
UI STUFF
================
*/

class TimeSeries {
    constructor(capacity, keys) {
        this.series = {};
        for (let key of keys) {
            this.series[key] = new Array(capacity).fill(0);
        }
        this.capacity = capacity;
        this.ix = 0;
    }
    set(key, value) {
        this.series[key][this.ix] = value;
    }
    get(key) {
        return this.series[key][this.ix];
    }
    advance() {
        this.ix += 1;
        this.ix = this.ix % this.capacity;
    }
    *iterator() {
        for (let key of Object.keys(this.series)) {
            yield { key, iterator: this.seriesIterator(key) };
        }
    }
    *seriesIterator(key) {
        for (let i = (this.ix + 1) % this.capacity; i != this.ix; i = (i + 1) % this.capacity) {
            yield this.series[key][i];
        }
    }
}

let vars = new TimeSeries(20, ['fps', 'E', 'L']);
let expfactor = 1;

function get_dt() {
    let this_time = Date.now();
    let dt = null;
    if (last_time != null) {
        dt = (this_time - last_time) / 1000;
    }
    else {
        dt = timestep / 1000;
    }
    last_time = this_time;
    return dt;
}

function main() {
    let dt = timestep / 1000;
    let canvas = document.getElementById('visuals');
    let timeseries = document.getElementById('timeseries');
    let frame = 0;
    setInterval(function() {
        draw(canvas, get_dt() / 1)
        let ctx = timeseries.getContext('2d');
        ctx.setTransform(1, 0, 0, 1, 0, 0);
        ctx.clearRect(0, 0, timeseries.width, timeseries.height);
        ctx.setTransform(1, 0, 0, -1, 0, timeseries.height/2);
        ctx.strokeStyle = "#c0c0c0";
        ctx.beginPath();
        ctx.moveTo(0, 0); ctx.lineTo(100, 0);
        ctx.stroke();
        for (let ts of vars.iterator()) {
            ctx.strokeStyle = ({fps: '#000000', E: '#ffff00', L: '#0000ff'})[ts.key];
            ctx.beginPath();
            let iter = ts.iterator;
            let v = iter.next().value;
            ctx.moveTo(0, 10 * Math.log(1+Math.abs(v)) * Math.sign(v));
            let x = 5;
            for (let v of iter) {
                ctx.lineTo(x, 10 * Math.log(1+Math.abs(v)) * Math.sign(v));
                x += 5
            }
            ctx.stroke();
        }
        console.log('E', vars.get('E'))
        frame++;
        if (frame % expfactor == 0) {
            vars.advance();
        }
    }, timestep);
}

function draw(canvas, dt) {
    if (dt > 0.1) {
        dt = 0.1;
    }
    let ctx = canvas.getContext('2d');
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    let scale = 0.25 * Math.min(canvas.width, canvas.height);
    ctx.setTransform(scale, 0, 0, -scale, canvas.width/2, canvas.height/2);
    ctx.lineWidth = 1 / scale;
    ctx.strokeStyle = '#000000';
    vars.set('fps', 1/dt);
    step(dt);
    render(ctx, dt);
}

function round(arr, places) {
    if (!places) {
        places = 2;
    }
    let pow = Math.pow(10, places);
    let out = [...arr];
    for (let i = 0; i < out.length; i++) {
        out[i] = Math.round(out[i]*pow)/pow;
    }
    return out;
}

function render(ctx, dt) {
    ctx.beginPath();
    ctx.moveTo(0, 0);
    let poss = pendulums(position);
    for (let i = 0; i < poss.length; i += 2) {
        ctx.lineTo(poss[i], poss[i+1]);
    }
    ctx.stroke();
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.fillText("pos: " + round(position), 0, 10);
    ctx.fillText("vel: " + round(velocity), 0, 20);
    ctx.fillText("acc: " + round(acceleration), 0, 30);
    ctx.fillText(iter, 0, 40);
    ctx.fillText("posxy: " + round(poss), 0, 50);
    let L_val = L(position, velocity);
    vars.set('L', L_val);
    ctx.fillText("L:  " + L_val, 0, 60);
    //ctx.fillText("L': " + L_closed(position, velocity), 0, 70);
    ctx.fillText("E: " + vars.get('E'), 0, 80);
    ctx.fillText("dt: " + dt, 0, 90);
    ctx.fillText("err: " + round([err], 5), 0, 100);
}

/*
================
CALCULUS
================
*/

function swap(arr, ix, val) {
    let copy = [...arr];
    copy.splice(ix, 1, val);
    return copy;
}

function differential(f, pos) {
    function df(x) {
        let prev = f(swap(x, pos, x[pos] - eps));
        let next = f(swap(x, pos, x[pos] + eps));
        for (let i = 0; i < prev.length; i++) {
            prev[i] = (next[i] - prev[i]) / (2*eps);
        }
        return prev;
    }
    return df;
}
function J(f) {
    let df = [];
    function j(x) {
        if (df.length == 0) {
            for (let i = 0; i < x.length; i++) {
                df.push(differential(f, i));
            }
        }
        let Jmat = [];
        for (let i = 0; i < x.length; i++) {
            Jmat.push(df[i](x));
        }
        return mattrans(Jmat);
    }
    return j;
}
function gradient(f) {
    function nabla(x) {
        let grad = [];
        for (let i = 0; i < x.length; i++) {
            let prev = f(swap(x, i, x[i] - eps))[0];
            let next = f(swap(x, i, x[i] + eps))[0];
            grad.push((next - prev) / (2*eps));
        }
        return grad;
    }
    return nabla;
}
function diff(f) {
    return function(x) { return (f(x+eps) - f(x-eps))/(2*eps) };
}

/*
================
PHYSICS
================
*/

let g = -9.82;
function transform_lagrangian(L, f) {
    function transformed(position, velocity) {
        let new_pos = f(position);
        let new_vel = [];
        for (let i = 0; i < new_pos.length; i++) {
            new_vel.push(0.0);
        }
        for (let i = 0; i < position.length; i++) {
            let dnew_pos = differential(f, i)(position);
            for (let j = 0; j < new_pos.length; j++) {
                new_vel[j] += velocity[i] * dnew_pos[j];
            }
        }
        return L(new_pos, new_vel);
    }
    return transformed;
}
function objects_langrangian(masses) {
    function L(poss, vels) {
        let out = 0;
        Energy = 0;
        for (let i = 0; i < masses.length; i++) {
            let m = masses[i];
            let px = poss[2*i]; let py = poss[2*i+1];
            let vx = vels[2*i]; let vy = vels[2*i+1];
            let K = 0.5 * m * (vx*vx + vy*vy);
            let P = -m * g * py;
            out += K - P;
            Energy += K + P;
        }
        console.log('Energy', Energy);
        vars.set('E', Energy);
        return [out];
    }
    return L;
}
function multi_pendulum(lengths) {
    function f(thetas) {
        let out = [];
        let x = 0;
        let y = 0;
        let theta = 0;
        for (let i = 0; i < lengths.length; i++) {
            theta = 0;
            let l = lengths[i]; theta += thetas[i];
            x += l * Math.sin(theta);
            y += l * Math.cos(theta);
            out.push(x); out.push(y);
        }
        return out;
    }
    return f;
}

/*
function L_closed(theta, dtheta) {
    l = 0.5;
    m = 1;
    l2 = Math.pow(l, 2);
    let K = 0.5 * (m+m) * l2 * Math.pow(dtheta[0], 2) +
            0.5 * m * l2 * Math.pow(dtheta[1], 2) +
            m * l2 * dtheta[0] * dtheta[1] * Math.cos(theta[0] - theta[1]);
    let P = -(m+m) * g * l * Math.cos(theta[0]) -
            m * g * l * Math.cos(theta[1]);
    let L = K - P;
    vars.set('E', K + P);
    console.log('KP', K, P)
    console.log('theta', theta, dtheta)
    console.log('L_inner', L)
    return [L];
}
*/

/*
================
SIMULATION
================
*/

pendulums = multi_pendulum([0.5, 0.5]);
let L = transform_lagrangian(objects_langrangian([1, 1]), pendulums);
let position = [Math.PI * 3/4, Math.PI*5/8];
let velocity = [0, 0];
let acceleration = [0, 0]; // To be computed
let first_acceleration = true;
//L = L_closed;
let iter = "N/A"; // ugly hack to print out # of newton iterations as debug value
let err = "N/A"; // ugly hack to print out accuracy to which the E-L equation was solved as a debug value
function step(dt) {
    for (let i = 0; i < position.length; i++) {
        position[i] += dt * velocity[i];
    }
    iter = 0;
    err = 1;
    while (err > 0.0001 && iter < 10) {
        err = 0;
        //console.log('L', L(position, velocity));
        // dL/dv = momentum
        let momentum = (p, v) => gradient((vel) => L(p, vel))(v);
        //console.log('momentum', momentum(position, velocity));
        // dL/dx = force
        let Force = (p, v) => gradient((pos) => L(pos, v))(p);
        //console.log('F', Force(position, velocity));
        // d/dt dL/dv = d/dt momentum = mass * acceleration
        let mass_acceleration = (acc) => differential((t) => {
            let vel_ = [...velocity];
            for (let i = 0; i < vel_.length; i++) vel_[i] += t[0] * acc[i];
            return momentum(position, vel_);
        }, 0)([0])
        //console.log('ma', mass_acceleration(acceleration));
        // Euler-Langrange: d/dt dL/dv = dL/dx    aka   ma = F
        let EulerLagrangeDifference = (acc) => {
            let F = Force(position, velocity);
            let ma = mass_acceleration(acc);
            let d = [];
            for (let i = 0; i < F.length; i++) {
                d.push(F[i] - ma[i]);
            }
            return d;
        }
        // solve using Newton's method
        let offset = EulerLagrangeDifference(acceleration);
        //console.log('offset', offset);
        let slope = J(EulerLagrangeDifference)(acceleration);
        //console.log(slope.toString())
        //console.log('slope', slope);
        let invSlope = matinv(slope);
        //let vec = offset;
        //let vecTest = matvecmul(slope, (matvecmul(invSlope, vec)));
        //for (let i = 0; i < vecTest.length; i++) {
        //    vecTest[i] -= vec[i];
        //}
        //console.log('vectest', vecTest);
        let shift = matvecmul(invSlope, offset);
        for (let i = 0; i < acceleration.length; i++) {
            acceleration[i] -= shift[i];
        }
        let error = EulerLagrangeDifference(acceleration);
        console.log(error);
        for (let i = 0; i < error.length; i++) {
            err += Math.abs(error[i]);
        }
        iter += 1;
    }
    first_acceleration = false;
    for (let i = 0; i < position.length; i++) {
        velocity[i] += dt * acceleration[i];
    }
}