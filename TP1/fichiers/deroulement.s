	.file	"deroulement.c"
	.section	.rodata
.LC1:
	.string	"sum %ld %.3f sec\n"
	.text
	.globl	main
	.type	main, @function
main:
.LFB2:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%rbx
	subq	$72, %rsp
	.cfi_offset 3, -24
	movl	$16384, %edi
	call	malloc
	movq	%rax, -48(%rbp)
	movl	$0, -28(%rbp)
	jmp	.L2
.L3:
	movl	-28(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-48(%rbp), %rax
	leaq	(%rdx,%rax), %rbx
	call	random
	movq	%rax, (%rbx)
	addl	$1, -28(%rbp)
.L2:
	cmpl	$2047, -28(%rbp)
	jle	.L3
	movl	$0, -36(%rbp)
	jmp	.L4
.L9:
	leaq	-64(%rbp), %rax
	movl	$0, %esi
	movq	%rax, %rdi
	call	gettimeofday
	movq	$0, -24(%rbp)
	movl	$0, -32(%rbp)
	jmp	.L5
.L8:
	movl	$0, -28(%rbp)
	jmp	.L6
.L7:
	movl	-28(%rbp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	movq	-48(%rbp), %rax
	addq	%rdx, %rax
	movq	(%rax), %rax
	addq	%rax, -24(%rbp)
	addl	$1, -28(%rbp)
.L6:
	cmpl	$2047, -28(%rbp)
	jle	.L7
	addl	$1, -32(%rbp)
.L5:
	cmpl	$1999999, -32(%rbp)
	jle	.L8
	leaq	-80(%rbp), %rax
	movl	$0, %esi
	movq	%rax, %rdi
	call	gettimeofday
	movq	-80(%rbp), %rdx
	movq	-64(%rbp), %rax
	subq	%rax, %rdx
	movq	%rdx, %rax
	imulq	$1000000, %rax, %rax
	movq	-72(%rbp), %rcx
	movq	-56(%rbp), %rdx
	subq	%rdx, %rcx
	movq	%rcx, %rdx
	addq	%rdx, %rax
	cvtsi2ssq	%rax, %xmm0
	movss	.LC0(%rip), %xmm1
	divss	%xmm1, %xmm0
	unpcklps	%xmm0, %xmm0
	cvtps2pd	%xmm0, %xmm0
	movq	-24(%rbp), %rax
	movq	%rax, %rsi
	movl	$.LC1, %edi
	movl	$1, %eax
	call	printf
	addl	$1, -36(%rbp)
.L4:
	cmpl	$9, -36(%rbp)
	jle	.L9
	movl	$0, %eax
	addq	$72, %rsp
	popq	%rbx
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE2:
	.size	main, .-main
	.section	.rodata
	.align 4
.LC0:
	.long	1232348160
	.ident	"GCC: (Debian 4.8.2-14) 4.8.2"
	.section	.note.GNU-stack,"",@progbits
