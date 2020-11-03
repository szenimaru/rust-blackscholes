//reference
//https://swdrsker.hatenablog.com/entry/2018/07/19/112450

fn normal(x: f64) -> f64 {
    let PI = 3.14159265358979323846264338327950288f64;
    let tmp1 = 2.0 * PI;
    let y = tmp1.sqrt();
    let z = x*x;
    let tmp2 = (-z)/2.0;
    let w = tmp2.exp();
    return w/y;
}

fn cdf(x: f64) -> f64 {
    let mut tmp_x = -100.0;
    let mut y = 0.0;
    let delta = 0.001;
    while tmp_x <= x {
        y += normal(tmp_x)*delta;
        tmp_x += delta;
    }
    return y;
}

fn BS_Call(S0: f64, sigma: f64,r: f64,T: f64, K: f64) -> f64{
    let tmp1 = S0 / K;
    let d1 = ( tmp1.ln() + ( r + sigma*sigma / 2.0 ) * T ) / ( sigma * T.sqrt() );
    let d2 = ( tmp1.ln() + ( r - sigma*sigma / 2.0 ) * T ) / ( sigma * T.sqrt() );
    let tmp2 = -r * T;
    let BS_C = S0 * cdf(d1) - K * tmp2.exp() * cdf(d2);
    return BS_C;
}

fn BS_Put(S0: f64, sigma: f64, r: f64, T: f64, K: f64) -> f64{
    let tmp1 = S0 / K;
    let d1 = ( tmp1.ln() + ( r + sigma*sigma / 2.0 ) * T ) / ( sigma * T.sqrt() );
    let d2 = ( tmp1.ln() + ( r - sigma*sigma / 2.0 ) * T ) / ( sigma * T.sqrt() );
    let tmp2 = -r * T;
    let BS_P = K * tmp2.exp() * cdf(-d2) - S0 * cdf(-d1);
    return BS_P;
}

fn main() {
    println!("Black Scholes");

    let S0 = 100.0;
    let sigma = 0.12;
    let r=0.005;
    let T=0.5;

    //println!("{} {}",BS_Call(S0,sigma,r,T,K),BS_Put(S0,sigma,r,T,K));

    let mut K = 50.0;
    while K < 200.0 {
        println!("{} {} {}",K,BS_Call(S0,sigma,r,T,K),BS_Put(S0,sigma,r,T,K));
        K += 1.0;
    }

    let tmp1 = -r * T;
    println!("{} {}",BS_Call(S0,sigma,r,T,K) + K * tmp1.exp(), BS_Put(S0,sigma,r,T,K)+S0);

}
