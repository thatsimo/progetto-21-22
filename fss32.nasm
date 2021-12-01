; ---------------------------------------------------------
; FSS con istruzioni SSE a 32 bit
; ---------------------------------------------------------

%include "sseutils32.nasm"
; max_vector_32(VECTOR x, int n, type* max)

section .data
section .bss
section .text

global max_vector_32
    x equ 8
	n equ 12
	max equ 16

max_vector_32:
	start
	mov eax,[ebp+x]
	mov edi,[ebp+n]
	sub edi,3
	movaps 	xmm0,[eax]
	mov 	esi, 4

fori_max:
    movaps 	xmm1,[eax+esi*4]
	minps 	xmm0,xmm1
	add	esi,4

	cmp 	esi,edi
	jl 	fori_max
	add 	edi,3


	movaps 	xmm1,xmm0
	shufps 	xmm1,xmm0,00001110b
	minps 	xmm0,xmm1
	movaps 	xmm1,xmm0
	shufps 	xmm1,xmm1,00000001b
	minps 	xmm0,xmm1

	cmp 	esi,edi
	jge 	end_max

forino_max:
	movss 	xmm1,[eax+esi*4]
	minss 	xmm0, xmm1
	add	esi,1
    cmp esi,edi
	jl	forino_max

end_max:
	mov 	eax,[ebp+max]
	movss 	[eax], xmm0

	stop
	
	
	
;euclidian_distance_32(MATRIX x, int offset, VECTOR y, int d,type* dist)
section .data
section .bss
section .text

global euclidian_distance_32
    x equ 8
    offset equ 12
    y equ 16
    d equ 20
    dist equ 24
euclidian_distance_32:

    start
    
    
    mov     eax,[ebp+x]
    mov     ebx,[ebp+offset]
    mov     ecx,[ebp+y]
    imul    ebx,4
    add     eax,ebx            ; indirizzo del vettore target
    
    mov     edi,[ebp+d]
    sub     edi,3
    
    xorps   xmm0,xmm0
    mov     esi,0
fori_euc:   
    movaps xmm1,[eax+esi*4]
    movaps xmm2,[ecx+esi*4]
    subps  xmm1,xmm2
    mulps  xmm1,xmm1
    addps  xmm0,xmm1
    
    add    esi,4
    cmp    esi,edi
    jl fori_euc
    
    add   edi,3
    
    haddps xmm0,xmm0
    haddps xmm0,xmm0
        
    cmp   esi,edi
    jge   end_euc
    
forino_euc:
    
    movss xmm1,[eax+esi*4]
    movss xmm2,[ecx+esi*4]
    subss xmm1,xmm2
    mulss xmm1,xmm1
    addss xmm0,xmm1
    add esi,1
    cmp esi,edi
    jl forino_euc
    
end_euc:    
    sqrtss xmm0,xmm0
    
    mov     eax,[ebp+dist]
    movss   [eax],xmm0
    
    stop

section .data
	
section .bss
section .text

global eval_f_32
       x1 equ 8
       d1 equ 12
       c1 equ 16
       offset1 equ 20
       quad equ 24
       scalar equ 28
       
eval_f_32:
       

    start

    xorps xmm0,xmm0     ; quad
    xorps xmm1,xmm1     ; scalar
    
    mov eax,[ebp+x1]
    mov ebx,[ebp+offset1]
    imul ebx,4
    add eax,ebx
    
    
    mov ebx,[ebp+c1]
    mov edi,[ebp+d1]
    sub edi,3
    mov esi,0
    
fori_exp:

    movaps xmm2,[eax+esi*4]
    movaps xmm3,xmm2
     
    mulps xmm2,xmm2
    addps xmm0,xmm2
    
    movaps xmm4,[ebx+esi*4]
    mulps xmm3,xmm4
    addps xmm1,xmm3
     
     
    add esi,4
    cmp esi,edi
    jl fori_exp
    
    
    add edi,3
    
    
    haddps xmm0,xmm0
    haddps xmm0,xmm0
    
    haddps xmm1,xmm1
    haddps xmm1,xmm1
    
    
    cmp esi,edi
    jge end_exp
    
forino_exp:
    
    movss xmm2,[eax+esi*4]
    movss xmm3,xmm2
     
    mulss xmm2,xmm2
    addss xmm0,xmm2
    
    movss xmm4,[ebx+esi*4]
    mulss xmm3,xmm4
    addss xmm1,xmm3
    
    inc esi
    cmp esi,edi
    jl forino_exp
    
    
 end_exp:
 
     mov eax,[ebp+quad]
     mov ebx,[ebp+scalar]
     
     movss [eax],xmm0
     movss [ebx],xmm1
     
     
     stop 
   
    

