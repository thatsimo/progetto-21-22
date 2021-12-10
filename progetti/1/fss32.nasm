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

	dim equ 4
	p equ 4
	UNROLL equ 4
max_vector_32:

	start
	mov eax,[ebp+x]
	mov edi,[ebp+n]
	sub edi,3


	movaps 	xmm0,[eax]
	mov 	esi, 4

fori:
    	movaps 	xmm1,[eax+esi*4]
	minps 	xmm0,xmm1

	add	esi,4
	cmp 	esi,edi
	jl 	fori


	add 	edi,3


	movaps 	xmm1,xmm0
	shufps 	xmm1,xmm0,00001110b
	minps 	xmm0,xmm1
	movaps 	xmm1,xmm0
	shufps 	xmm1,xmm1,00000001b
	minps 	xmm0,xmm1

	cmp 	esi,edi
	jge 	end

forino:
	movss 	xmm1,[eax+esi*4]
	minss 	xmm0, xmm1
	add	esi,1
	jmp 	forino

end:
	mov 	eax,[ebp+max]
	movss 	[eax], xmm0

	stop
	
	
	
;euclidian_distance_32(MATRIX x, int offset, VECTOR y, int d,type* dist)
%include "sseutils32.nasm"
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
fori:   
    movaps xmm1,[eax+esi*4]
    movaps xmm2,[ecx+esi*4]
    subps  xmm1,xmm2
    mulps  xmm1,xmm1
    addps  xmm0,xmm1
    
    add    esi,4
    cmp    esi,edi
    jl fori
    
    add   edi,3
    
    haddps xmm0,xmm0
    haddps xmm0,xmm0
        
    cmp   esi,edi
    jge   end
    
forino:
    
    movss xmm1,[eax+esi*4]
    movss xmm2,[ecx+esi*4]
    subss xmm1,xmm2
    mulss xmm1,xmm1
    addss xmm0,xmm1
    add esi,1
    cmp esi,edi
    jl forino
    
end:    
    sqrtss xmm0
    
    mov     eax,[ebp+dist]
    movss   [eax],xmm0
    
    stop


