%include "sseutils64.nasm"
; min_vector_64(VECTOR x, int n, type* max)

section .data
section .bss
section .text

global min_vector_64
	; rdi = x
    ; rsi = n
    ; rdx = max
    UNROLL_MAX equ 8 
min_vector_64 :
	start

    vmovupd ymm0,[rdi]
    mov     rax,0
    sub     rsi,UNROLL_MAX-1

fori_min:
    vmovupd ymm1,[rdi+rax*8]
    vminpd  ymm0,ymm1

    vmovupd ymm1,[rdi+rax*8+32]
    vminpd  ymm0,ymm1

    add     rax,UNROLL_MAX
    cmp     rax,rsi
    jl      fori_min

    add     rsi,UNROLL_MAX-1

    vperm2f128 ymm1,ymm0,ymm0,00010001b
    vminpd     ymm0,ymm1
    vmovupd    xmm1,xmm0
    vshufpd    xmm0,xmm0,01b
    vminpd     xmm0,xmm1
    
    cmp    rax,rsi
    jge    end_min

forino_min:

    vmovq   xmm1,[rdi+rax*8]
    vminsd  xmm0,xmm1
    inc     rax
    cmp     rax,rsi
    jl      forino_min

end_min:

    vmovq   [rdx],xmm0

		

	stop
	

;vector_sum_64(MATRIX x, int offset, int n,VECTOR v);

section .data
section .bss
section .text

global vector_sum_64
    ; rdi = x
    ; rsi = offset
    ; rdx = n
    ; rcx = v

    UNROLL_VS equ 8 

	; controllare sopra 

vector_sum_64:
    start

    imul rsi,8
    add rdi,rsi

    sub rdx,UNROLL_VS-1
    mov rax,0


fori_ws:
    vmovupd ymm0,[rdi+rax*8]
    vaddpd   ymm0,[rcx+rax*8]
    vmovupd [rdi+rax*8],ymm0

    vmovupd ymm0,[rdi+rax*8+32]
    vaddpd   ymm0,[rcx+rax*8+32]
    vmovupd [rdi+rax*8+32],ymm0

    add rax,UNROLL_VS
    cmp rax,rdx
    jl fori_ws

    add rdx,UNROLL_VS-1

    cmp rax,rdx
    jge end_vs

forino_ws:

    vmovq xmm0,[rdi+rax*8]
    vaddsd xmm0,[rcx+rax*8]
    vmovupd [rdi+rax*8],xmm0

    inc rax
    cmp rax,rdx
    jl forino_ws

end_vs:

    stop


;euclidian_distance_64(MATRIX x, int offset, VECTOR y, int d,type* dist)
section .data
section .bss
section .text

global euclidian_distance_64
    ; rdi = x
    ; rsi = offset
    ; rdx = y
    ; rcx = d
    ; r8 = dist

    UNROLL_EUC equ 8

euclidian_distance_64:
    start
    
    imul    rsi,8
    add     rdi,rsi

    mov     rax,0
    sub     rcx,UNROLL_EUC-1

    vxorpd ymm0,ymm0

fori_euc:
    vmovupd ymm1,[rdi+rax*8]
    vmovupd ymm2,[rdx+rax*8]
    vsubpd  ymm1,ymm2
    vmulpd  ymm1,ymm1
    vaddpd  ymm0,ymm1

    vmovupd ymm1,[rdi+rax*8+32]
    vmovupd ymm2,[rdx+rax*8+32]
    vsubpd  ymm1,ymm2
    vmulpd  ymm1,ymm1
    vaddpd  ymm0,ymm1
    
    add rax,UNROLL_EUC
    cmp rax,rcx
    jl fori_euc


    add rcx,UNROLL_EUC-1

    vhaddpd ymm0,ymm0,ymm0
    vperm2f128 ymm2,ymm0,ymm0,00010001b
    vaddpd xmm0,xmm2

    cmp rax,rcx
    jge end_euc

forino_euc:
    vmovq   xmm1,[rdi+rax*8]
    vmovq   xmm2,[rdx+rax*8]
    vsubsd  xmm1,xmm2
    vmulsd  xmm1,xmm1
    vaddsd  xmm0,xmm1
    inc     rax
    cmp     rax,rcx
    jl      forino_euc

end_euc:
    vsqrtsd   xmm0,xmm0
    vmovq   [r8],xmm0    
    stop


;eval_f_64(VECTOR x,int d,VECTIOR c, int offset,type* quad,type* scalar)

section .data
section .bss
section .text
global eval_f_64
    ; rdi = x
    ; rsi = d
    ; rdx = c
    ; rcx = offset
    ; r8 = quad
    ; r9 = scalar
    UNROLL_EVALF equ 8

eval_f_64:

    start
    
    imul rcx,8
    add rdi,rcx

    mov rax,0

    vxorpd ymm0,ymm0
    vxorpd ymm1,ymm1

    sub rsi,UNROLL_EVALF-1


fori_f:
    vmovupd ymm4,[rdi+rax*8]
    vmovupd ymm2,ymm4
    vmulpd  ymm2,ymm2
    vaddpd  ymm0,ymm2

    vmovupd ymm3,[rdx+rax*8]
    vmulpd  ymm3,ymm4
    vaddpd  ymm1,ymm3


    vmovupd ymm4,[rdi+rax*8+32]
    vmovupd ymm2,ymm4
    vmulpd  ymm2,ymm2
    vaddpd  ymm0,ymm2

    vmovupd ymm3,[rdx+rax*8+32]
    vmulpd  ymm3,ymm4
    vaddpd  ymm1,ymm3

    add rax,UNROLL_EVALF
    cmp rax,rsi
    jl fori_f

    add rsi,UNROLL_EVALF-1

    vhaddpd ymm0,ymm0,ymm0
    vperm2f128 ymm2,ymm0,ymm0,00010001b
    vaddpd xmm0,xmm2
    

    vhaddpd ymm1,ymm1,ymm1
    vperm2f128 ymm2,ymm1,ymm1,00010001b
    vaddpd xmm1,xmm2

    cmp rax,rsi
    jge end_f

forino_f:

    vmovq   xmm2,[rdi+rax*8]
    vmovq   xmm3,xmm2
    vmulsd  xmm3,xmm3
    vaddsd  xmm0,xmm3

    vmovq   xmm3,[rdx+rax*8]
    vmulsd  xmm2,xmm3
    vaddsd  xmm1,xmm2

    inc rax
    cmp rax,rsi
    jl forino_f
    
end_f:
    vmovq [r8],xmm0
    vmovq [r9],xmm1

    stop
    



;compute_avg_64(MATRIX x, int np, int d, VECTOR c,type den, VECTOR ris)
section .data
section .bss
section .text

global compute_avg_64
    ; rdi = x 
    ; rsi = np
    ; rdx = d
    ; rcx = c
    ; r8 = ris
    ; xmm0 = den
    UNROLL_WA equ 8
compute_avg_64:

    start
    
     
    vbroadcastsd    ymm0,xmm0
       
    mov rax,0                       ; i=0

fori_wa:
    mov          rbx,0              ; j=0
    vbroadcastsd ymm1,[rcx+rax*8]
    vdivpd       ymm1,ymm0
    mov          r10,rdx
    imul         r10,rax            ; i*d
    sub          rdx,UNROLL_WA-1
forj_wa:
    mov         r11,r10
    add         r11,rbx             ; i*d+j
    vmovupd     ymm2,[rdi+r11*8]
    vmulpd      ymm2,ymm1
    vaddpd      ymm2,[r8+rbx*8]
    vmovupd     [r8+rbx*8],ymm2

    vmovupd     ymm2,[rdi+r11*8+32]
    vmulpd      ymm2,ymm1
    vaddpd      ymm2,[r8+rbx*8+32]
    vmovupd     [r8+rbx*8+32],ymm2

    add         rbx,UNROLL_WA
    cmp         rbx,rdx
    jl          forj_wa

    add         rdx,UNROLL_WA-1

    cmp         rbx,rdx
    jge         end_forj_wa

forino_wa:
    mov         r11,r10 
    add         r11,rbx             ; i*d+j
    vmovq       xmm2,[rdi+r11*8]
    vmulsd      xmm2,xmm1
    vaddsd      xmm2,[r8+rbx*8]
    vmovq       [r8+rbx*8],xmm2

    inc         rbx
    cmp         rbx,rdx
    jl          forino_wa

end_forj_wa:

    inc         rax
    cmp         rax,rsi
    jl          fori_wa

    stop
