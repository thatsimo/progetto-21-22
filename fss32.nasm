%include "sseutils32.nasm"

; min_vector_32(VECTOR x, int n, type* min)
section .data
section .bss
section .text

global min_vector_32
    x equ 8
	n equ 12
	min equ 16
    UNROLL_MIN equ 8

min_vector_32:
	start
	mov     eax,[ebp+x]             ; x
	mov     edi,[ebp+n]             ; n
	sub     edi,UNROLL_MIN-1        ; gestione vettore non multiplo
	movups 	xmm0,[eax]              ; primi quattro elementi
	mov 	esi, 4                  ; i=4

fori_min:
    movups 	xmm1,[eax+esi*4]        ; xmm1<-x[...]
	minps 	xmm0,xmm1               ; confronto xmm0 xmm1

    movups 	xmm1,[eax+esi*4+16]     ; UNROLL    
	minps 	xmm0,xmm1

	add	    esi,UNROLL_MIN          ; i+=UNROLL

	cmp 	esi,edi                 ; i<n-7?
	jl 	    fori_min

	add 	edi,UNROLL_MIN-1        ; ripristino n

	movups 	xmm1,xmm0               ; riduzione vettore
	shufps 	xmm1,xmm0,00001110b
	minps 	xmm0,xmm1
	movups 	xmm1,xmm0
	shufps 	xmm1,xmm1,00000001b
	minps 	xmm0,xmm1               ; primo elemento di xmm0
                                    ; per la parte scalare
	cmp 	esi,edi                 ; i<n? 
	jge 	end_min                 ; gestione caso lunghezza multipla di 8

forino_min:                         ; replica codice per versione scalare
	movss 	xmm1,[eax+esi*4]        
	minss 	xmm0, xmm1
	add	    esi,1
    cmp     esi,edi
	jl	    forino_min

end_min:                            ; caricamento del minimo in memoria
	mov 	eax,[ebp+min]
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
    UNROLL_EUC equ 8

euclidian_distance_32:
    start    
    mov     eax,[ebp+x]         ; x
    mov     ebx,[ebp+offset]    ; offset
    mov     ecx,[ebp+y]         ; y
    imul    ebx,4               ; offset in versione byte
    
    add     eax,ebx             ; indirizzo del vettore target nella matrice
    
    mov     edi,[ebp+d]         ; d
    sub     edi,UNROLL_EUC-1    ; gestione vettore non multiplo di 8
    xorps   xmm0,xmm0           ; ret=0
    mov     esi,0               ; i=0

fori_euc:   
    movups xmm1,[eax+esi*4]     ; xmm1<- x[...]
    movups xmm2,[ecx+esi*4]     ; xmm2<- x[...]
    subps  xmm1,xmm2            ; (x-y)
    mulps  xmm1,xmm1            ; (x-y)^2
    addps  xmm0,xmm1            ; ret+=(x-y)^2

    movups xmm1,[eax+esi*4+16]  ; UNROLL
    movups xmm2,[ecx+esi*4+16]
    subps  xmm1,xmm2
    mulps  xmm1,xmm1
    addps  xmm0,xmm1
    
    add    esi,UNROLL_EUC       ; i+=8

    cmp    esi,edi              ; i<n-7?
    jl     fori_euc
    
    add    edi,UNROLL_EUC        ; i+=8
    
    haddps xmm0,xmm0            ; riduzione di xmm0
    haddps xmm0,xmm0
    cmp   esi,edi               ; caso vettore multiplo di 8
    jge   end_euc
    
forino_euc:                     ; replica scalare della sezione precedente
    movss xmm1,[eax+esi*4]
    movss xmm2,[ecx+esi*4]
    subss xmm1,xmm2
    mulss xmm1,xmm1
    addss xmm0,xmm1
    add esi,1
    
    cmp esi,edi
    jl forino_euc
    
end_euc:                        
    sqrtss xmm0,xmm0            ; sqrt(ret)
    mov     eax,[ebp+dist]      ; caricamento in memoria 
    movss   [eax],xmm0
    stop
    
;eval_f_32(MATRIX x, int d, VECTOR c, int offset,type* quad, type* scalar);
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
       UNROLL_F equ 8
       
eval_f_32:
    start

    xorps xmm0,xmm0         ; quad=0
    xorps xmm1,xmm1         ; scalar=0
    
    mov eax,[ebp+x1]        ; x
    mov ebx,[ebp+offset1]   ; offset
    imul ebx,4              ; offset in versione byte
    add eax,ebx             ; eax <- indirizzo vettore target
    
    mov ebx,[ebp+c1]        ; c
    mov edi,[ebp+d1]        ; d
    sub edi,UNROLL_F-1      ; gestione vettore non multiplo di 8
    mov esi,0               ; i=0
    
fori_exp:
    movups  xmm2,[eax+esi*4]    ; xmm2 <- x[...]
    movups  xmm3,xmm2           ; duplicazione xmm2 in xmm3
     
    mulps   xmm2,xmm2           ; x^2
    addps   xmm0,xmm2           ; quad+=x^2    
    
    movups  xmm4,[ebx+esi*4]    ; xmm4 <- c[...]
    mulps   xmm3,xmm4           ; x*c 
    addps   xmm1,xmm3           ; scalar+=x*c

    movups  xmm2,[eax+esi*4+16] ; UNROLL
    movups  xmm3,xmm2
     
    mulps   xmm2,xmm2
    addps   xmm0,xmm2
    
    movups  xmm4,[ebx+esi*4+16]
    mulps   xmm3,xmm4
    addps   xmm1,xmm3
     
    add     esi,UNROLL_F
    
    cmp     esi,edi              ; i<n-7?
    jl      fori_exp
    
    add edi,UNROLL_F-1          ; ripristino lunghezza 
    
    haddps xmm0,xmm0            ; riduzione quad
    haddps xmm0,xmm0
    
    haddps xmm1,xmm1            ; riduzione scalar
    haddps xmm1,xmm1
    
    cmp esi,edi                 ; i<n?
    jge end_exp
    
forino_exp:                     ; versione scalare della sezione precedente
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
    
 end_exp:                       ; caricamento in memoria
    mov eax,[ebp+quad]
    mov ebx,[ebp+scalar]
     
    movss [eax],xmm0
    movss [ebx],xmm1
    stop 

;compute_avg_32(MATRIX x, int np, int d, VECTOR c,type den, VECTOR ris)
section .data
section .bss
section .text

global compute_avg_32
    x_w equ 8
    np_w equ 12
    d_w equ 16
    c_w equ  20
    den_w equ 24
    ris_w equ 28
   
compute_avg_32:
    start
    
    mov eax,[ebp+x_w]              ; x
    mov ebx,[ebp+c_w]              ; c
    mov ecx,[ebp+d_w]              ; d
    movd xmm5,esp
    mov esp,[ebp+ris_w]
    
    movss xmm0,[ebp+den_w]
    shufps xmm0,xmm0,00000000b    ; den

    mov esi,0                     ; i=0
    
foriw: 
    mov edi,0                     ; j=0
    movss xmm7,[ebx+esi*4]

    shufps xmm7,xmm7,00000000b
    divps xmm7,xmm0
    mov edx,1

    imul edx,ecx
    imul edx,esi

    sub ecx,7

forjw:
    add edx,edi
    movups xmm1,[eax+edx*4]

    sub edx,edi
    mulps xmm1,xmm7

    addps xmm1,[esp+edi*4]
    movups [esp+edi*4],xmm1

    add edx,edi
    movups xmm1,[eax+edx*4+16]

    sub edx,edi
    mulps xmm1,xmm7

    addps xmm1,[esp+edi*4+16]
    movups [esp+edi*4+16],xmm1
    
    add edi, 8

    cmp edi,ecx
    jl forjw
    
    add ecx,7
    
    cmp edi,ecx
    jge endforjw
    
forinow:
    add edx,edi
    movss xmm1,[eax+edx*4]

    sub edx,edi
    mulss xmm1,xmm7

    addss xmm1,[esp+edi*4]
    movups [esp+edi*4],xmm1
    
    inc edi

    cmp edi,ecx
    jl forinow
    
endforjw:
    inc esi

    cmp esi,[ebp+np_w]
    jl foriw
    movd esp,xmm5   
    stop 

;vector_sum32(MATRIX x, int offset, int n,VECTOR v);
section .data
section .bss
section .text

global vector_sum_32
	x_vs equ 8
	offset_vs equ 12
	d_vs equ 16
	v_vs equ 20

vector_sum_32:
	start

	mov esi,0			; i=0
	mov eax,[ebp+x_vs]		; x
	mov ecx,[ebp+offset]		; offset

	imul ecx,4			
	
    add eax,ecx			; x[i*d]
	mov ebx,[ebp+v_vs]		; v
	mov edi,[ebp+d_vs]		; d
	sub edi,3			; d-3
	
fori_vs:
	movups xmm0,[eax+esi*4]
	movups xmm1,[ebx+esi*4]
	
	addps  xmm0,xmm1
	movups [eax+esi*4],xmm0
	
	add esi,4
	
    cmp esi,edi
	jl fori_vs

	add edi,3

	cmp esi,edi
	jge end_vs

forino_vs:
	movss xmm0,[eax+esi*4]
	movss xmm1,[ebx+esi*4]
	
	addss xmm0,xmm1
	
    movss [eax+esi*4],xmm0
	
    inc esi
	
    cmp esi,edi
	jl forino_vs

end_vs:
	stop