#[cfg(test)]
mod test {
    use pasta_curves::group::ff::PrimeField;
    use poseidon_pasta_curve::{Fr, Poseidon};

    #[test]
    fn test2val() {
        let mut poseidon = Poseidon::new(2);
        poseidon.inputs[0] = Fr::from_str_vartime("11").unwrap();
        poseidon.inputs[1] = Fr::from_str_vartime("22").unwrap();
        poseidon.run();
        assert!(poseidon.out.eq(&Fr::from_str_vartime(
            "9534981536951832176089863541273014817120680825980532463463767768652007267989"
        )
        .unwrap()));
        println!("{:?}", poseidon.out);
    }
}
