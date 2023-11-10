# calc_cor_of_cor works

    Code
      calc_cor_of_cor(testing_ex_list)
    Message
      No parallel processing has been detected
      Computing intra-study correlations using 328 genes
      Converting correlations to Z scores...
    Output
                GSE39582      READ      COAD
      GSE39582 1.0000000 0.4550780 0.5031691
      READ     0.4550780 1.0000000 0.7991573
      COAD     0.5031691 0.7991573 1.0000000

# compare_networks works

    Code
      compare_networks(testing_memb_list[[1]], testing_memb_list[[2]])
    Output
          adjRand cor_of_cor  cor_coph overlap
      1 0.3584286  0.3293828 0.1466221      19

