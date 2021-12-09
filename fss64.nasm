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

    vmovapd ymm0,[rdi]
    mov     rax,0
    sub     rsi,3

fori_min:
    vmovapd ymm1,[rdi+rax*8]
    vminpd  ymm0,ymm1
    add     rax,4
    cmp     rax,rsi
    jl      fori_min

    add     rsi,3

    vperm2f128 ymm1,ymm0,ymm0,00010001b
    vminpd     ymm0,ymm1
    vmovapd    xmm1,xmm0
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

    sub rdx,3
    mov rax,0


fori_ws:
    vmovapd ymm0,[rdi+rax*8]
    vaddpd   ymm0,[rcx+rax*8]
    vmovapd [rdi+rax*8],ymm0

    add rax,4
    cmp rax,rdx
    jl fori_ws

    add rdx,3

    cmp rax,rdx
    jge end_vs

forino_ws:

    vmovq xmm0,[rdi+rax*8]
    vaddsd xmm0,[rcx+rax*8]
    vmovapd [rdi+rax*8],xmm0

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
    sub     rcx,3

    vxorpd ymm0,ymm0

fori_euc:
    vmovapd ymm1,[rdi+rax*8]
    vmovapd ymm2,[rdx+rax*8]
    vsubpd  ymm1,ymm2
    vmulpd  ymm1,ymm1
    vaddpd  ymm0,ymm1
    
    add rax,4
    cmp rax,rcx
    jl fori_euc


    add rcx,3

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
    jl      fori_euc

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

    sub rsi,3


fori_f:
    vmovapd ymm4,[rdi+rax*8]
    vmovapd ymm2,ymm4
    vmulpd  ymm2,ymm2
    vaddpd  ymm0,ymm2

    vmovapd ymm3,[rdx+rax*8]
    vmulpd  ymm3,ymm4
    vaddpd  ymm1,ymm3

    add rax,4
    cmp rax,rsi
    jl fori_f

    add rsi,3

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
    ; r8 = den
    ; r9 = ris
    UNROLL_W equ 8
compute_avg_64:

    start
    mov rax,rdi     ; x
    mov rbx,rcx     ; c
    mov r10,rdx     ; d
    mov r11,rsi     ; np
    

    vbroadcastsd ymm0,[r8]
    
    mov rsi,0
    
    
foriw:
    
    mov             rdi,0
    vbroadcastsd    ymm7,[rbx+rsi*8]
    
    
    vdivpd  ymm7,ymm0
    
    mov     rdx,r10
    imul    rdx,rsi

    sub     r10,UNROLL_W-1
forjw:
    
    add     rdx,rdi
    vmovupd ymm1,[rax+rdx*8]
    sub     rdx,rdi
    vmulpd  ymm1,ymm7
    vaddpd  ymm1,[r9+rdi*8]
    vmovupd [r9+rdi*8],ymm1


    add     rdx,rdi
    vmovupd ymm1,[rax+rdx*8+32]
    sub     rdx,rdi
    vmulpd  ymm1,ymm7
    vaddpd  ymm1,[r9+rdi*8+32]
    vmovupd [r9+rdi*8+32],ymm1
    
    add     rdi,UNROLL_W
    cmp     rdi,r10
    jl      forjw

    add     r10,UNROLL_W-1

    cmp     rdi,r10
    jge     endforjw
    
forinow:
    add     rdx,rdi
    vmovq   xmm1,[rax+rdx*8]
    sub     rdx,rdi
    vmulsd  xmm1,xmm7
    vaddsd  xmm1,[r9+rdi*8]
    vmovq   [r9+rdi*8],xmm1
    
    inc rdi
    cmp rdi,r10
    jl forinow
    

endforjw:

    inc rsi 
    cmp rsi,r11
    jl  foriw


    stop
