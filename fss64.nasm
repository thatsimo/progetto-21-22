%include "sseutils64.nasm"

; min_vector_64(VECTOR x, int n, type* max)

section .data
section .bss
section .text

global min_vector_64
	x equ 8
	n equ 12
	max equ 16
    UNROLL_MAX equ 4 ; ripristinare a 8
min_vector_64 :
	start
	mov     rax,[rbp+x]             ; x
	mov     rdi,[rbp+n]             ; n
	sub     rdi,UNROLL_MAX-1        ; gestione vettore non multiplo
	vmovapd ymm0,[rax]              ; primi quattro elementi (in teoria)
	mov 	rsi, 4                  ; i=4

fori_min:
    vmovapd ymm1,[rax+rsi*8]        ; xmm1<-x[...]
	vminpd 	ymm0,ymm1               ; confronto ymm0 ymm1

   ; vmovapd ymm1,[rax+rsi*8+16]     ; UNROLL    
;	vminpd 	ymm0,ymm1

	add	    rsi,UNROLL_MAX          ; i+=UNROLL

	cmp 	rsi,rdi                 ; i<n-7?
	jl 	    fori_min


	add 	rdi,UNROLL_MAX-1        ; ripristino n


	vmovapd ymm1,ymm0               ; riduzione vettore
	vshufpd ymm1,ymm0,00001110b
	vminpd 	ymm0,ymm1
	vmovapd ymm1,ymm0
	vshufpd ymm1,ymm1,00000001b
	vminpd 	ymm0,ymm1               


	cmp 	rsi,rdi                 ; i<n? 
	jge 	end_min                 ; gestione caso lunghezza multipla di 8

forino_min:                         ; replica codice per versione scalare double
	vmovq 	xmm1,[rax+rsi*8]    	; muovo un double in xmm1. Nota : la v sembrerebbe non necessaria ma dovrebbe impedire il context switch evitando di incorrere nella penalità    
	vminpd 	xmm0, xmm1				; la v è presente per lo stesso ragionamento di sopra
	add	    rsi,1
    cmp     rsi,rdi
	jl	    forino_min

end_min:                            ; caricamento del minimo in memoria
	mov 	rax,[rbp+max]
	vmovq 	[rax], xmm0			

	stop
	


;vector_sum64(MATRIX x, int offset, int n,VECTOR v);

section .data
section .bss
section .text

global vector_sum_64
    x_vs equ 8
    offset_vs equ 12
    n_vs equ 16
    v_vs equ 20

    UNROLL_VS equ 4 ; ripristinare a 8 in seguito

	; controllare sopra 

vector_sum_64:
    start

    mov     rax,[ebp+x_vs]          ; x
    mov     rbx,[ebp+offset_vs]     ; offset
    imul    rbx,4                   ; porta offset a versione byte 
    add     rax,rbx                 ; porta l'indice alla posizione del vettore target

    mov     rbx,[rbp+v_vs]          ; v

    mov     rdi,[rbp+n_vs]          ; n
    sub     rdi,UNROLL_VS-1         ; unroll

    mov     rsi,0                   ; i=0
fori_vs:

    vmovapd  ymm0,[rbx+rsi*8]        ; v[...]
    vaddpd   ymm0,[rax+rsi*8]        ; somma v[...] con x[...]
    vmovapd  [rax+rsi*8],ymm0        ; carica il risultato su x[...] in memoria

   ; vmovapd  ymm0,[rbx+rsi*8+16]     ; UNROLL
   ; vaddpd   ymm0,[rax+rsi*8+16]
   ; vmovapd  [rax+rsi*8+16],ymm0

   ; add     rsi,UNROLL_VS           ; i+=8

    cmp     rsi,rdi                 ; i<n-7?
    jl      fori_vs


    add     rdi,UNROLL_VS-1         ; ripristino n

    cmp     rsi,rdi                 ; i<n?
    jge     end_vs
    
forino_vs:                          ; gestione caso vettore non multiplo di 8
    vmovq   xmm0,[rbx+rsi*8]  
    vmovq   xmm1, [rax+rsi*8] 		; vedere se si possa evitare lo spostamento nel registro e farlo direttamente dalla memoria  
    vaddpd  xmm0,xmm1 
    vmovq   [rcx+rsi*8],xmm0

    inc     rsi
    cmp     rsi,rdi
    jl      forino_vs

end_vs:
    stop


