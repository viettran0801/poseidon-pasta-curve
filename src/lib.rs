mod constants;

use crate::constants::{poseidon_c, poseidon_m, poseidon_p, poseidon_s};
use num::{BigInt, Num};
use pasta_curves::group::ff::PrimeField;

pub type Fr = pasta_curves::pallas::Scalar;

#[derive(Clone, Debug)]
struct Sigma {
    pub inn: Fr,
    pub out: Fr,
}

impl Sigma {
    pub fn default() -> Sigma {
        Sigma::new()
    }

    pub fn new() -> Sigma {
        Sigma {
            inn: Fr::zero(),
            out: Fr::zero(),
        }
    }

    pub fn run(&mut self) {
        let i2 = self.inn.mul(&self.inn);
        let i4 = i2.mul(&i2);
        self.out = i4.mul(&self.inn);
    }
}

#[derive(Clone, Debug)]
struct Ark {
    pub inn: Vec<Fr>,
    pub out: Vec<Fr>,
    pub t: usize,
    pub c: Vec<Fr>,
    pub r: usize,
}

impl Ark {
    pub fn default() -> Ark {
        Ark {
            inn: vec![],
            out: vec![],
            t: 0,
            c: vec![],
            r: 0,
        }
    }

    pub fn new(t: usize, c: Vec<Fr>, r: usize) -> Ark {
        Ark {
            inn: vec![Fr::zero(); t],
            out: vec![Fr::zero(); t],
            t,
            c,
            r,
        }
    }

    pub fn run(&mut self) {
        for i in 0..self.t {
            self.out[i] = self.inn[i].add(&self.c[i + self.r]);
        }
    }
}

#[derive(Clone, Debug)]
struct Mix {
    pub inn: Vec<Fr>,
    pub out: Vec<Fr>,
    pub t: usize,
    pub m: Vec<Vec<Fr>>,
}

impl Mix {
    pub fn default() -> Mix {
        Mix {
            inn: vec![],
            out: vec![],
            t: 0,
            m: vec![],
        }
    }

    pub fn new(t: usize, m: Vec<Vec<Fr>>) -> Mix {
        Mix {
            inn: vec![Fr::zero(); t],
            out: vec![Fr::zero(); t],
            t,
            m,
        }
    }

    pub fn run(&mut self) {
        for i in 0..self.t {
            let mut lc = Fr::zero();
            for j in 0..self.t {
                lc = lc.add(&self.m[j][i].mul(&self.inn[j]));
            }
            self.out[i] = lc.clone();
        }
    }
}

#[derive(Clone, Debug)]
struct MixLast {
    pub inn: Vec<Fr>,
    pub out: Fr,
    pub t: usize,
    pub m: Vec<Vec<Fr>>,
    pub s: usize,
}

impl MixLast {
    pub fn default() -> MixLast {
        MixLast {
            inn: vec![],
            out: Fr::zero(),
            t: 0,
            m: vec![],
            s: 0,
        }
    }

    pub fn new(t: usize, m: Vec<Vec<Fr>>, s: usize) -> MixLast {
        MixLast {
            inn: vec![Fr::zero(); t],
            out: Fr::zero(),
            t,
            m,
            s,
        }
    }

    pub fn run(&mut self) {
        let mut lc = Fr::zero();
        for j in 0..self.t {
            lc = lc.add(&self.m[j][self.s].mul(&self.inn[j]));
        }
        self.out = lc.clone()
    }
}

#[derive(Clone, Debug)]
struct MixS {
    pub inn: Vec<Fr>,
    pub out: Vec<Fr>,
    pub t: usize,
    pub s: Vec<Fr>,
    pub r: usize,
}

impl MixS {
    pub fn default() -> MixS {
        MixS {
            inn: vec![],
            out: vec![],
            t: 0,
            s: vec![],
            r: 0,
        }
    }

    pub fn new(t: usize, s: Vec<Fr>, r: usize) -> MixS {
        MixS {
            inn: vec![Fr::zero(); t],
            out: vec![Fr::zero(); t],
            t,
            s,
            r,
        }
    }

    pub fn run(&mut self) {
        let mut lc = Fr::zero();
        for i in 0..self.t {
            lc = lc.add(&self.s[(self.t * 2 - 1) * self.r + i].mul(&self.inn[i]));
        }

        self.out[0] = lc;

        for i in 1..self.t {
            let tmp = self.inn[0].mul(&self.s[(self.t * 2 - 1) * self.r + self.t + i - 1]);
            self.out[i] = tmp.add(&self.inn[i]);
        }
    }
}

#[derive(Clone, Debug)]
struct PoseidonEx {
    pub inputs: Vec<Fr>,
    pub initial_state: Fr,
    pub out: Vec<Fr>,
    pub n_inputs: usize,
    pub n_outs: usize,
}

impl PoseidonEx {
    pub fn new(n_inputs: usize, n_outs: usize) -> PoseidonEx {
        PoseidonEx {
            inputs: vec![Fr::zero(); n_inputs],
            initial_state: Fr::zero(),
            out: vec![Fr::zero(); n_outs],
            n_inputs,
            n_outs,
        }
    }

    pub fn run(&mut self) {
        let n_round_p: Vec<usize> = vec![
            56, 57, 56, 60, 60, 63, 64, 63, 60, 66, 60, 65, 70, 60, 64, 68,
        ];
        let t = self.n_inputs + 1;
        let n_round_f: usize = 8;
        let n_round_p = n_round_p[t - 2];

        let c: Vec<Fr> = poseidon_c(t)
            .iter()
            .map(|x| {
                Fr::from_str_vartime(
                    BigInt::from_str_radix(&x[2..], 16)
                        .unwrap()
                        .to_string()
                        .as_str(),
                )
                .unwrap()
            })
            .collect();

        let s: Vec<Fr> = poseidon_s(t)
            .iter()
            .map(|x| {
                Fr::from_str_vartime(
                    BigInt::from_str_radix(&x[2..], 16)
                        .unwrap()
                        .to_string()
                        .as_str(),
                )
                .unwrap()
            })
            .collect();

        let m: Vec<Vec<Fr>> = poseidon_m(t)
            .iter()
            .map(|x| {
                x.iter()
                    .map(|xx| {
                        Fr::from_str_vartime(
                            BigInt::from_str_radix(&xx[2..], 16)
                                .unwrap()
                                .to_string()
                                .as_str(),
                        )
                        .unwrap()
                    })
                    .collect::<Vec<Fr>>()
            })
            .collect();

        let p: Vec<Vec<Fr>> = poseidon_p(t)
            .iter()
            .map(|x| {
                x.iter()
                    .map(|xx| {
                        Fr::from_str_vartime(
                            BigInt::from_str_radix(&xx[2..], 16)
                                .unwrap()
                                .to_string()
                                .as_str(),
                        )
                        .unwrap()
                    })
                    .collect::<Vec<Fr>>()
            })
            .collect();

        let mut ark = vec![Ark::default(); n_round_f];
        let mut sigma_f = vec![vec![Sigma::default(); t]; n_round_f];
        let mut sigma_p = vec![Sigma::default(); n_round_p];
        let mut mix = vec![Mix::default(); n_round_f - 1];
        let mut mix_s = vec![MixS::default(); n_round_p];
        let mut mix_last = vec![MixLast::default(); self.n_outs];

        ark[0] = Ark::new(t, c.clone(), 0);

        for j in 0..t {
            if j > 0 {
                ark[0].inn[j] = self.inputs[j - 1].clone();
            } else {
                ark[0].inn[j] = self.initial_state.clone();
            }
        }

        ark[0].run();

        for r in 0..(n_round_f / 2 - 1) {
            for j in 0..t {
                sigma_f[r][j] = Sigma::new();
                if r == 0 {
                    sigma_f[r][j].inn = ark[0].out[j].clone();
                } else {
                    sigma_f[r][j].inn = mix[r - 1].out[j].clone();
                }
                sigma_f[r][j].run();
            }

            ark[r + 1] = Ark::new(t, c.clone(), (r + 1) * t);
            for j in 0..t {
                ark[r + 1].inn[j] = sigma_f[r][j].out.clone();
            }

            ark[r + 1].run();

            mix[r] = Mix::new(t, m.clone());
            for j in 0..t {
                mix[r].inn[j] = ark[r + 1].out[j].clone();
            }
            mix[r].run();
        }

        for j in 0..t {
            sigma_f[n_round_f / 2 - 1][j] = Sigma::new();
            sigma_f[n_round_f / 2 - 1][j].inn = mix[n_round_f / 2 - 2].out[j].clone();
            sigma_f[n_round_f / 2 - 1][j].run();
        }

        ark[n_round_f / 2] = Ark::new(t, c.clone(), (n_round_f / 2) * t);
        for j in 0..t {
            ark[n_round_f / 2].inn[j] = sigma_f[n_round_f / 2 - 1][j].out.clone();
        }
        ark[n_round_f / 2].run();

        mix[n_round_f / 2 - 1] = Mix::new(t, p.clone());
        for j in 0..t {
            mix[n_round_f / 2 - 1].inn[j] = ark[n_round_f / 2].out[j].clone();
        }

        mix[n_round_f / 2 - 1].run();

        for r in 0..n_round_p {
            sigma_p[r] = Sigma::new();
            if r == 0 {
                sigma_p[r].inn = mix[n_round_f / 2 - 1].out[0].clone();
            } else {
                sigma_p[r].inn = mix_s[r - 1].out[0].clone();
            }

            sigma_p[r].run();

            mix_s[r] = MixS::new(t, s.clone(), r);

            for j in 0..t {
                if j == 0 {
                    mix_s[r].inn[j] = sigma_p[r].out.add(&c[(n_round_f / 2 + 1) * t + r]);
                } else {
                    if r == 0 {
                        mix_s[r].inn[j] = mix[n_round_f / 2 - 1].out[j].clone();
                    } else {
                        mix_s[r].inn[j] = mix_s[r - 1].out[j].clone();
                    }
                }
                mix_s[r].run();
            }
        }

        for r in 0..(n_round_f / 2 - 1) {
            for j in 0..t {
                sigma_f[n_round_f / 2 + r][j] = Sigma::new();
                if r == 0 {
                    sigma_f[n_round_f / 2 + r][j].inn = mix_s[n_round_p - 1].out[j].clone();
                } else {
                    sigma_f[n_round_f / 2 + r][j].inn = mix[n_round_f / 2 + r - 1].out[j].clone();
                }
                sigma_f[n_round_f / 2 + r][j].run();
            }

            ark[n_round_f / 2 + r + 1] =
                Ark::new(t, c.clone(), (n_round_f / 2 + 1) * t + n_round_p + r * t);

            for j in 0..t {
                ark[n_round_f / 2 + r + 1].inn[j] = sigma_f[n_round_f / 2 + r][j].out.clone();
            }
            ark[n_round_f / 2 + r + 1].run();

            mix[n_round_f / 2 + r] = Mix::new(t, m.clone());

            for j in 0..t {
                mix[n_round_f / 2 + r].inn[j] = ark[n_round_f / 2 + r + 1].out[j].clone();
            }
            mix[n_round_f / 2 + r].run();
        }

        for j in 0..t {
            sigma_f[n_round_f - 1][j] = Sigma::new();
            sigma_f[n_round_f - 1][j].inn = mix[n_round_f - 2].out[j].clone();
            sigma_f[n_round_f - 1][j].run();
        }

        for i in 0..self.n_outs {
            mix_last[i] = MixLast::new(t, m.clone(), i);
            for j in 0..t {
                mix_last[i].inn[j] = sigma_f[n_round_f - 1][j].out.clone();
            }
            mix_last[i].run();
            self.out[i] = mix_last[i].out.clone();
        }
    }
}

pub struct Poseidon {
    pub inputs: Vec<Fr>,
    pub out: Fr,
    pub n_inputs: usize,
}

impl Poseidon {
    pub fn new(n_inputs: usize) -> Poseidon {
        Poseidon {
            inputs: vec![Fr::zero(); n_inputs],
            out: Fr::zero(),
            n_inputs,
        }
    }

    pub fn run(&mut self) {
        let mut p_ex = PoseidonEx::new(self.n_inputs, 1);
        p_ex.initial_state = Fr::zero();
        for i in 0..self.n_inputs {
            p_ex.inputs[i] = self.inputs[i];
        }
        p_ex.run();
        self.out = p_ex.out[0].clone();
    }
}
