#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.scripts_dir = "$projectDir"  // You can override this with --scripts_dir option

process DO_A {
    input:
    path script

    output:
    path '*', emit: a_output

    script:
    """
    bash $script
    """
}

process DO_B {
    input:
    path script

    output:
    path '*', emit: b_output

    script:
    """
    bash $script
    """
}

process DO_C {
    input:
    path script
    val a_done
    val b_done

    output:
    path '*', emit: c_output

    script:
    """
    bash $script
    """
}

workflow {
    println "Scripts directory: ${params.scripts_dir}"
    
    do_a_script = file("${params.scripts_dir}/do_a.sh")
    do_b_script = file("${params.scripts_dir}/do_b.sh")
    do_c_script = file("${params.scripts_dir}/do_c.sh")

    println "Script files:"
    println "A: ${do_a_script.name} (exists: ${do_a_script.exists()})"
    println "B: ${do_b_script.name} (exists: ${do_b_script.exists()})"
    println "C: ${do_c_script.name} (exists: ${do_c_script.exists()})"

    a_result = DO_A(do_a_script)
    b_result = DO_B(do_b_script)
    DO_C(do_c_script, a_result.a_output, b_result.b_output)
}