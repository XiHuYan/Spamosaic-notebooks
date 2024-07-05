## Mouse embryonic dataset
CUDA_VISIBLE_DEVICES=2 python run.py --exp e3 --task misar_E15-13-18 --epoch_num 2000
CUDA_VISIBLE_DEVICES=2 python run.py --task misar_E15-13-18 --act predict_all_latent_bc --init_model sp_00001999 --exp e1
CUDA_VISIBLE_DEVICES=2 python run.py --exp e1 --task misar_E15_18_18_modalMatch --epoch_num 2000
CUDA_VISIBLE_DEVICES=2 python run.py --task misar_E15_18_18_modalMatch --act predict_all_latent_bc --init_model sp_00001999 --exp e1

## Mouse postnatal brain RNA+ATAC
CUDA_VISIBLE_DEVICES=1 python run.py --exp e1 --task MB_RNA_ATAC --epoch_num 2000
CUDA_VISIBLE_DEVICES=1 python run.py --task MB_RNA_ATAC --act predict_all_latent_bc --init_model sp_00001999 --exp e1

## Mouse postnatal brain RNA+H3K4me3
CUDA_VISIBLE_DEVICES=1 python run.py --exp e1 --task MB_RNA_H3K4me3 --epoch_num 2000
CUDA_VISIBLE_DEVICES=1 python run.py --task MB_RNA_H3K4me3 --act predict_all_latent_bc --init_model sp_00001999 --exp e1

## Mouse postnatal brain RNA+H3K27me3
CUDA_VISIBLE_DEVICES=1 python run.py --exp e1 --task MB_RNA_H3K27me3 --epoch_num 2000
CUDA_VISIBLE_DEVICES=1 python run.py --task MB_RNA_H3K27me3 --act predict_all_latent_bc --init_model sp_00001999 --exp e1

## Mouse postnatal brain RNA+H3K27ac
CUDA_VISIBLE_DEVICES=1 python run.py --exp e1 --task MB_RNA_H3K27ac --epoch_num 2000
CUDA_VISIBLE_DEVICES=1 python run.py --task MB_RNA_H3K27ac --act predict_all_latent_bc --init_model sp_00001999 --exp e1
CUDA_VISIBLE_DEVICES=2 python run.py --exp e1 --task MB_RNA_H3K27ac_modalMatch --epoch_num 2000
CUDA_VISIBLE_DEVICES=2 python run.py --task MB_RNA_H3K27ac_modalMatch --act predict_all_latent_bc --init_model sp_00001999 --exp e1

## Mouse embryo
CUDA_VISIBLE_DEVICES=1 python run.py --exp e1 --task ME_RNA_ATAC --epoch_num 2000
CUDA_VISIBLE_DEVICES=1 python run.py --task ME_RNA_ATAC --act predict_all_latent_bc --init_model sp_00001999 --exp e1

## Lymph node
CUDA_VISIBLE_DEVICES=2 python run.py --exp e1 --task Lymph --epoch_num 2000
CUDA_VISIBLE_DEVICES=2 python run.py --task Lymph --act predict_all_latent_bc --init_model sp_00001999 --exp e1
CUDA_VISIBLE_DEVICES=2 python run.py --exp e1 --task Lymph_modalMatch --epoch_num 2000
CUDA_VISIBLE_DEVICES=2 python run.py --task Lymph_modalMatch --act predict_all_latent_bc --init_model sp_00001999 --exp e1

## Tonsil
CUDA_VISIBLE_DEVICES=2 python run.py --exp e1 --task Tonsil --epoch_num 2000
CUDA_VISIBLE_DEVICES=2 python run.py --task Tonsil --act predict_all_latent_bc --init_model sp_00001999 --exp e1
CUDA_VISIBLE_DEVICES=2 python run.py --exp e1 --task Tonsil_modalMatch --epoch_num 2000
CUDA_VISIBLE_DEVICES=2 python run.py --task Tonsil_modalMatch --act predict_all_latent_bc --init_model sp_00001999 --exp e1