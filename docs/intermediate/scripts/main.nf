#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process DO_A {
    output:
    path '*', emit: a_output

    script:
    """
    bash do_a.sh
    """
}

process DO_B {
    output:
    path '*', emit: b_output

    script:
    """
    bash do_b.sh
    """
}

process DO_C {
    input:
    val a_done
    val b_done

    output:
    path '*', emit: c_output

    script:
    """
    bash do_c.sh
    """
}

workflow {
    a_result = DO_A()
    b_result = DO_B()
    DO_C(a_result.a_output, b_result.b_output)
}